"""Numerical epicycle extraction by frequency analysis (NAFF), derivation Part 5.

Method (Laskar 1990, Icarus 88, `Laskar1990`; Laskar 1993, Physica D 67,
`Laskar1993`, in research/refs.bib).  For a uniformly sampled series f(t_n),
n = 0..N-1, with the weighted inner product <f,g> = sum_n w_n f_n conj(g_n),
where w is the Hann^p window  w_n ~ sin^(2p)(pi n/(N-1))  (Laskar's chi_p;
p = 1 is the plain Hann window), normalised to sum(w) = 1:

1. Locate the strongest peak of the windowed (zero-padded x4) FFT of the
   current residual.
2. Refine the frequency by maximising, with a bounded Brent search within
   +-1 FFT bin of the peak, the power of the projection of the residual on the
   candidate basis function *orthogonalised against all previously accepted
   basis functions* (Gram-Schmidt / QR in the weighted space).  This is
   Laskar's |<f, e^{i w t}>| maximisation with the orthogonalisation made
   explicit; for real data the candidate is the pair {cos wt, sin wt}, which
   removes the interference between the +w and -w lines of a real signal.
3. Subtract the projection (orthogonal subtraction) and repeat for K terms.
4. Final amplitudes: joint weighted least squares over all accepted
   frequencies (the same subspace the sequential Gram-Schmidt amplitudes
   span, but without error accumulation).

Window power: leakage of far lines falls off faster with larger p, but the
main lobe widens (half-width ~ (p+1) bins).  Measured on a 200-yr Mars
geocentric series (5-day sampling, 30 frequencies; scratch benchmark, also
asserted loosely in tests/test_naff.py): worst relative frequency / radius
error over the 10 largest analytic terms is 2e-10 / 5e-7 for p=3 and
7e-11 / 1e-7 for p=4 (p=2 is ~3e-9 / 5e-4 or worse on the weak terms).
Default p = 4 (Laskar's sin^8 window).  Note the sample-weighted fit is
least accurate at the record ends where the window vanishes: the pointwise
reconstruction error of the chain is ~5e-10 AU (max over the record) in that
benchmark; on the Eros 100-yr series it is 7.6e-7 AU at the worst (record ends)
but 3.1e-8 AU in the middle 80%, so report full-window and interior errors.

Frequency refinement: scipy's bounded Brent has tolerance sqrt(eps)*|x|, i.e.
~1e-8 relative in an absolute frequency, so the search is done on the offset
from the FFT peak, in nested windows (see :func:`_maximize`).  Without this
the Earth line was off by 3e-9 relative and the least-squares step produced
spurious nearby companion lines.  The optional ``sweeps`` backfitting passes
(re-maximising each frequency with all other terms projected out, aimed at
greedy-order bias; default 0 = plain Laskar NAFF) were found *not* to help in
benchmarks (they can let two lines collapse onto one another, and are slow),
so they are off by default.

Resolution limit (measured in tests/test_naff.py, see there for numbers):
two lines are resolved only if separated by more than ~2*pi/T (T = record
length); with the Hann^p window the practical limit is roughly (p+1) bins
(~(p+1)*2pi/T) before interference biases the frequencies.  Closer pairs are
returned as one (merged, biased) line plus, at most, a poorly conditioned
partner with large cancelling amplitudes; outputs stay finite.  Isolated
well-separated lines are recovered to ~1e-10 relative frequency accuracy
(measured 1.5e-10 .. 4.6e-10) on 100-period records (see test), limited by the flatness of the objective near
its maximum (sqrt(machine eps) of a bin).

Public API: :func:`naff_complex`, :func:`naff_vector`.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import fft, next_fast_len
from scipy.optimize import minimize_scalar

from .expansion import Epicycle, ellipse_to_circles


def hann_window(N: int, p: int = 4) -> np.ndarray:
    """Laskar's Hann^p window sin^(2p)(pi n/(N-1)), normalised to sum 1."""
    if N < 2:
        raise ValueError("need at least 2 samples")
    w = np.sin(np.pi * np.arange(N) / (N - 1)) ** (2 * p)
    return w / w.sum()


def _uniform_dt(t):
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.size < 8:
        raise ValueError("t must be a 1-D array with at least 8 samples")
    d = np.diff(t)
    dt = (t[-1] - t[0]) / (t.size - 1)
    if dt <= 0 or np.max(np.abs(d - dt)) > 1e-6 * dt:
        raise ValueError("t must be uniformly increasing")
    return dt


def _maximize(negpow, lo, hi, scale):
    """Minimise negpow on [lo, hi] by nested bounded Brent searches.

    scipy's bounded Brent has a relative tolerance sqrt(eps)*|x|, which for an
    absolute frequency x is far coarser than needed; so the search is done on
    the *offset* from the interval centre, then repeated in successively
    narrower windows (+-1e-2 and +-1e-4 bins; one bin is
    ``scale``) centred on the previous optimum.
    """
    c, h = 0.5 * (lo + hi), 0.5 * (hi - lo)
    best = None
    for half in (h, 1e-2 * scale, 1e-4 * scale):
        a, b = max(lo, c - half), min(hi, c + half)
        r = minimize_scalar(lambda d: negpow(c + d), bounds=(a - c, b - c), method="bounded",
                            options={"xatol": 1e-13 * scale})
        if best is None or r.fun <= best[1]:
            best = (c + r.x, r.fun)
        c = best[0]
    return best


def _engine(R, t, K, real, window_p, constant, pad, search_bins, rtol, cond_tol, sweeps=0):
    """Shared greedy projection pursuit.  R: (N, c) data, real or complex.

    Returns (freqs, coeffs, const, t_ref): coeffs shape (m, 2, c) for real
    (cos, sin coefficients) or (m, 1, c) for complex; const (c,) or None.
    """
    N, c = R.shape
    dt = _uniform_dt(t)
    t_ref = float(t[0])
    tau = np.asarray(t, dtype=float) - t_ref
    sw = np.sqrt(hann_window(N, window_p))[:, None]
    Rs = sw * R
    dtype = float if (real and not np.iscomplexobj(R)) else complex
    Rs = Rs.astype(dtype)
    norm0 = np.linalg.norm(Rs)
    Q = np.zeros((N, 0), dtype=dtype)
    cols = []  # scaled basis columns of accepted terms
    if constant:
        q0 = sw[:, 0].astype(dtype)
        cols_const = q0.copy()
        q0 = q0 / np.linalg.norm(q0)
        Q = q0[:, None]
        Rs = Rs - Q @ (Q.conj().T @ Rs)
    else:
        cols_const = None
    Np = next_fast_len(int(pad) * N)
    fgrid = 2 * np.pi * np.fft.fftfreq(Np, dt)
    bin_w = 2 * np.pi / (N * dt)
    nyq = np.pi / dt

    def basis(om):
        if real:
            return sw * np.stack([np.cos(om * tau), np.sin(om * tau)], axis=1)
        return sw * np.exp(1j * om * tau)[:, None]

    def projected(om, Q=None):
        Q = _Q[0] if Q is None else Q
        B = basis(om)
        s0 = np.linalg.norm(B, axis=0).max()
        if Q.shape[1]:
            B = B - Q @ (Q.conj().T @ B)
        U, s, _ = np.linalg.svd(B, full_matrices=False)
        keep = s > cond_tol * s0
        return U[:, keep], s0

    _Q = [Q]
    freqs = []
    for _ in range(K):
        if Q.shape[1] >= N - 2 or np.linalg.norm(Rs) <= rtol * norm0:
            break
        spec = np.sum(np.abs(fft(Rs, n=Np, axis=0)) ** 2, axis=1)
        if real:
            ok = (fgrid > 0)
        else:
            ok = np.ones_like(fgrid, bool)
        idx = np.flatnonzero(ok)
        om0 = fgrid[idx[np.argmax(spec[idx])]]
        lo, hi = om0 - search_bins * bin_w, om0 + search_bins * bin_w
        if real:
            lo, hi = max(lo, 1e-3 * bin_w), min(hi, nyq)
        else:
            lo, hi = max(lo, -nyq), min(hi, nyq)

        def negpow(om):
            U, _ = projected(om)
            if U.shape[1] == 0:
                return 0.0
            return -float(np.linalg.norm(U.conj().T @ Rs) ** 2)

        om, fbest = _maximize(negpow, lo, hi, bin_w)
        om = float(om)
        U, _ = projected(om)
        if U.shape[1] == 0 or -fbest <= (rtol * norm0) ** 2:
            break
        freqs.append(om)
        cols.append(basis(om))
        Q = np.concatenate([Q, U], axis=1)
        _Q[0] = Q
        Rs = Rs - U @ (U.conj().T @ Rs)

    # Backfitting sweeps: the greedy pass picks each line while weaker lines are
    # still unsubtracted (their window leakage biases it); re-maximise each
    # frequency against the data with all *other* terms projected out.
    Rfull = (sw * R).astype(dtype)
    for _ in range(sweeps if len(freqs) > 1 else 0):
        for j in range(len(freqs)):
            others = ([cols_const[:, None]] if constant else []) + [
                cols[i] for i in range(len(cols)) if i != j]
            Qo, _r = np.linalg.qr(np.concatenate(others, axis=1))
            Rj = Rfull - Qo @ (Qo.conj().T @ Rfull)

            def negpow_j(om):
                U, _ = projected(om, Qo)
                return -float(np.linalg.norm(U.conj().T @ Rj) ** 2) if U.shape[1] else 0.0

            lo, hi = freqs[j] - 0.5 * bin_w, freqs[j] + 0.5 * bin_w
            if real:
                lo, hi = max(lo, 1e-3 * bin_w), min(hi, nyq)
            else:
                lo, hi = max(lo, -nyq), min(hi, nyq)
            omn, fn = _maximize(negpow_j, lo, hi, bin_w)
            if fn <= negpow_j(freqs[j]):
                freqs[j] = float(omn)
                cols[j] = basis(freqs[j])

    m = len(freqs)
    if m == 0 and not constant:
        return [], np.zeros((0, 2 if real else 1, c), dtype), None, t_ref
    parts = ([cols_const[:, None]] if constant else []) + cols
    G = np.concatenate(parts, axis=1)
    coef, *_ = np.linalg.lstsq(G, (sw * R).astype(dtype), rcond=None)
    const = None
    if constant:
        const = coef[0]
        coef = coef[1:]
    per = 2 if real else 1
    return freqs, coef.reshape(m, per, c), const, t_ref


def naff_complex(z, t, K, window_p: int = 4, pad: int = 4, search_bins: float = 1.0,
                 rtol: float = 1e-13, cond_tol: float = 1e-6, sweeps: int = 0):
    """NAFF of a uniformly sampled complex series z(t) ~ sum_k A_k exp(i w_k (t - t[0])).

    z: (N,) complex (real input is accepted and yields +-w pairs); t: (N,)
    uniform times; K: max number of terms (fewer returned if the residual
    reaches rtol * ||z||_w or the basis is exhausted).  cond_tol: candidate
    basis functions whose Gram-Schmidt remainder is below cond_tol of their
    norm (indistinguishable from earlier ones) are rejected.  sweeps: number
    of optional backfitting passes re-optimising every frequency with the
    others projected out (aimed at greedy-order bias, but not found helpful and
    can merge lines; default 0 = plain Laskar NAFF, see module docstring).
    Returns list of (omega [rad per unit of t], complex amplitude A at t[0]) in
    order of extraction (decreasing power); omega may be negative.
    """
    z = np.asarray(z)
    if z.ndim != 1 or z.shape != np.shape(t):
        raise ValueError("z and t must be 1-D arrays of equal length")
    freqs, coef, _, _ = _engine(z[:, None].astype(complex), t, K, False, window_p, False,
                                pad, search_bins, rtol, cond_tol, sweeps)
    return [(w, complex(coef[j, 0, 0])) for j, w in enumerate(freqs)]


def naff_vector(xyz, t, K, window_p: int = 4, pad: int = 4, search_bins: float = 1.0,
                rtol: float = 1e-13, cond_tol: float = 1e-6, sweeps: int = 0, source: str = "naff",
                t_ref: float | None = None) -> list[Epicycle]:
    """NAFF of a real 3-D trajectory xyz (N,3) sampled at uniform times t (JD).

    The constant (mean) vector is fitted first (an omega = 0 fixed arm), then
    up to K frequencies w > 0 are found on the summed power of all three
    components; the final amplitudes come from one joint weighted least-squares
    fit, giving each frequency a complex 3-vector A with term
    Re(A exp(i w (t - t_ref))), converted to a prograde + retrograde circle
    pair by :func:`expansion.ellipse_to_circles`.  t_ref defaults to t[0]
    (phases are re-referenced if another t_ref is given).  Result is sorted by
    decreasing radius and evaluable by :func:`expansion.evaluate_chain`.

    Convention: each line w > 0 becomes the circle pair from
    :func:`expansion.ellipse_to_circles`: a prograde circle (omega = +w) and a
    retrograde circle (omega = -w), both with the ellipse-plane normal; the
    constant term is an omega = 0 arm.  Note a circle (omega, normal) traces
    the same curve as (-omega, -normal), so when comparing with analytic chains
    canonicalise first (e.g. flip to omega > 0, negating normal and swapping
    u/v, or compare via the sign of omega*normal).
    """
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or xyz.shape[0] != np.size(t):
        raise ValueError("xyz must have shape (N, 3) matching t")
    freqs, coef, const, t0 = _engine(xyz, t, K, True, window_p, True, pad,
                                     search_bins, rtol, cond_tol, sweeps)
    if t_ref is None:
        t_ref = t0
    out: list[Epicycle] = []
    r = float(np.linalg.norm(const))
    if r > 0:
        u = const / r
        helper = np.eye(3)[np.argmin(np.abs(u))]
        v = np.cross(np.cross(u, helper), u)
        v /= np.linalg.norm(v)
        out.append(Epicycle(r, 0.0, 0.0, np.cross(u, v), u, v, t_ref, source, 0))
    for j, w in enumerate(freqs):
        a, b = coef[j, 0], coef[j, 1]       # a cos + b sin = Re((a - i b) e^{i w tau})
        A = (a - 1j * b) * np.exp(1j * w * (t_ref - t0))
        out.extend(ellipse_to_circles(A, w, t_ref, source, j + 1))
    out.sort(key=lambda ep: -ep.radius)
    return out
