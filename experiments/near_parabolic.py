"""Near-parabolic minimal-N study (review item 1) and NAFF line table for Keplerian Mars.

Run: uv run python experiments/near_parabolic.py   (about 2-3 min)
Writes experiments/out/near_parabolic.json and experiments/out/naff_lines_mars.json and
prints a LaTeX-ready table of each.

Part A. For z/a = sum_k c_k exp(ikM) (perifocal, focus at origin) compare, per eccentricity e
and tolerance tol (units of a, position error = |z - z_N|/a, complex modulus):
  N_cert   : smallest N with the paper's Eq.(16) bound epicycles.bounds.tail_bound(e, N) <= tol
             (None for e = 1: the bound is infinite there).
  N_actual : smallest N with sup_M |z - z_N|/a <= tol (and stays below for all larger N).
Method for N_actual.  Because z(M=0) = 1 - e exactly (periapsis) and c_0 = -3e/2, the tail at
M = 0 is  T0(N) = (1 + e/2) - sum_{k=1..N}(c_k + c_-k)  (no infinite sum needed).  If every
c_{+-k} with k > N is >= 0 then every term of the tail is maximal at M = 0 and
T0(N) = sum_{|k|>N}|c_k| = sup_M |tail|, i.e. the coefficient tail sum IS the sup.  Signs are
checked numerically for all k <= Kmax (reported as "signs") and the sup is cross-checked by a
dense sampled max over M (log-spaced near M = 0 plus uniform) against the exact Kepler solution.
Where T0 has not reached tol by Kmax the remaining N is extrapolated with the limiting law
T0 ~ 3 C N^(-2/3) (flagged "extrapolated").
Part B. Tail exponent for e = 1: c_{+-k} = J'_k(k)/k ~ C k^(-5/3), C = 3^(1/6) Gamma(2/3)/(2^(1/3) pi).
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import gamma, jv

from epicycles.bodies import EARTH, MARS
from epicycles.bounds import rho, tail_bound
from epicycles.expansion import geocentric_chain
from epicycles.frames import geocentric
from epicycles.kepler import solve_kepler
from epicycles.naff import naff_vector

OUT = Path(__file__).parent / "out"
ES = [0.5, 0.89, 0.99, 0.999, 0.9999, 0.999999, 1.0]
TOLS = [1e-2, 1e-6]
KMAX = 10_000_000
CHUNK = 1_000_000
C_LIM = 3 ** (1 / 6) * gamma(2 / 3) / (2 ** (1 / 3) * np.pi)   # J'_k(k) ~ C_LIM k^(-2/3)


def coeff_pair(e, k):
    """(c_k, c_-k) for integer array k >= 1; e = 1 is the limit J'_k(k)/k for both."""
    k = np.asarray(k, dtype=float)
    if e == 1.0:
        c = (jv(k, k) - jv(k + 1, k)) / k
        return c, c.copy()
    s = np.sqrt(1.0 - e * e)
    jk, jk1 = jv(k, k * e), jv(k + 1, k * e)
    return ((1.0 + s) / e * jk - jk1) / k, (e / (1.0 + s) * jk - jk1) / k


def n_cert(e, tol):
    """Smallest N with tail_bound(e, N) <= tol (monotone in N; exponential + bisection)."""
    if e >= 1.0:
        return None
    hi = 1
    while tail_bound(e, hi) > tol:
        hi *= 2
    lo = hi // 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if tail_bound(e, mid) > tol:
            lo = mid
        else:
            hi = mid
    return hi


def stream_tail(e, tols, kmax=KMAX, grid=True):
    """Stream coefficients; return first-crossing N per tol, sign info, tail checkpoints."""
    total = 1.0 + e / 2.0                     # sum_{k != 0} c_k = z(0) - c_0
    S, k0 = 0.0, 1
    cross = {t: None for t in tols}
    last_above = {t: 0 for t in tols}
    minc, firstneg = np.inf, None
    ck, tk = [], []                           # checkpoints (N, T0) on a ~geometric grid
    nxt = 1.0
    stop_below = min(tols) / 10
    while k0 <= kmax:
        k = np.arange(k0, min(k0 + CHUNK, kmax + 1), dtype=float)
        cp, cm = coeff_pair(e, k)
        s = cp + cm
        T = total - (S + np.cumsum(s))
        m = min(cp.min(), cm.min())
        minc = min(minc, m)
        if firstneg is None and m < 0:
            firstneg = int(k[np.argmax((cp < 0) | (cm < 0))])
        for t in tols:
            above = np.flatnonzero(np.abs(T) > t)
            if len(above):
                last_above[t] = int(k[above[-1]])
        if grid:
            idx = np.unique(np.minimum(np.searchsorted(k, nxt * 1.15 ** np.arange(0, 60)), len(k) - 1))
            for i in idx:
                if k[i] >= nxt:
                    ck.append(k[i]); tk.append(T[i]); nxt = k[i] * 1.15
        S += s.sum()
        k0 = int(k[-1]) + 1
        if abs(T[-1]) < stop_below:
            break
    for t in tols:
        cross[t] = last_above[t] + 1 if (abs(T[-1]) <= t) else None
    return cross, minc, firstneg, int(k[-1]), float(T[-1]), np.array(ck), np.array(tk)


def exact_z(e, M):
    if e < 1.0:
        E = solve_kepler(M, e)
        return (np.cos(E) - e) + 1j * np.sqrt(1 - e * e) * np.sin(E)
    lo, hi = np.zeros_like(M), np.full_like(M, 2 * np.pi)   # M in [0, pi]: E - sin E = M monotone
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        big = mid - np.sin(mid) > M
        hi = np.where(big, mid, hi); lo = np.where(big, lo, mid)
    return np.cos(0.5 * (lo + hi)) - 1.0 + 0j


def sampled_tail_max(e, N, nM=4000):
    """max over an M grid of |z - z_N| (complex modulus), plus the value at M = 0 and argmax M."""
    M = np.unique(np.concatenate([np.geomspace(1e-9, np.pi, 1500), np.linspace(0, np.pi, nM)]))
    cp, cm = coeff_pair(e, np.arange(1, N + 1))
    c0 = -1.5 * e
    zN = np.full(M.shape, c0, dtype=complex)
    ks = np.arange(1, N + 1)
    for j in range(0, len(M), 200):
        Mj = M[j:j + 200]
        ang = np.outer(Mj, ks)
        zN[j:j + 200] += np.cos(ang) @ (cp + cm) + 1j * (np.sin(ang) @ (cp - cm))
    err = np.abs(exact_z(e, M) - zN)
    i = int(np.argmax(err))
    return float(err[i]), float(M[i]), float(err[0])


def part_a():
    rows = []
    for e in ES:
        t0 = time.time()
        cross, minc, firstneg, kdone, Tlast, ck, tk = stream_tail(e, TOLS)
        row = {"e": e, "rho": None if e >= 1 else rho(e), "Kmax_done": kdone,
               "T0_at_Kmax": Tlast, "min_coeff_all_k": float(minc),
               "all_coeffs_nonneg": bool(minc >= 0), "first_negative_k": firstneg,
               "N_cert": {}, "N_actual": {}, "extrapolated": {}, "sampled": {}}
        for tol in TOLS:
            key = f"{tol:g}"
            row["N_cert"][key] = n_cert(e, tol)
            N = cross[tol]
            ex = N is None
            if ex:   # extrapolate with T0 ~ T0(Kmax) (N/Kmax)^(-2/3)
                N = int(round(kdone * (Tlast / tol) ** 1.5))
            row["N_actual"][key] = int(N)
            row["extrapolated"][key] = bool(ex)
            if not ex and N <= 30000 and e < 1.0 or (e == 1.0 and not ex and N <= 30000):
                mx, Mmx, T_0 = sampled_tail_max(e, N)
                mx1, _, _ = sampled_tail_max(e, N - 1) if N > 1 else (None, None, None)
                row["sampled"][key] = {"N": N, "sampled_max": mx, "argmax_M": Mmx, "tail_at_M0": T_0,
                                       "sampled_max_at_N_minus_1": mx1}
        rows.append(row)
        print(f"e={e}: done in {time.time()-t0:.0f}s  Kmax={kdone}  nonneg={row['all_coeffs_nonneg']}", flush=True)
    return rows


def part_b():
    """Tail exponent / constant for e = 1 from the streamed coefficients."""
    _, _, _, _, _, ck, tk = stream_tail(1.0, [1e-12], kmax=KMAX)
    ok = (ck >= 1e3) & (ck <= 1e7)
    p_tail = np.polyfit(np.log(ck[ok]), np.log(tk[ok]), 1)
    k = np.unique(np.round(np.geomspace(1e3, 1e7, 40))).astype(int)
    cp, _ = coeff_pair(1.0, k)
    p_c = np.polyfit(np.log(k), np.log(cp), 1)
    kk = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
    cc, _ = coeff_pair(1.0, kk)
    return {"C_closed_form": float(C_LIM), "coef_fit_exponent": float(p_c[0]),
            "coef_fit_C": float(np.exp(p_c[1])), "coef_exponent_expected": -5 / 3,
            "tail_fit_exponent": float(p_tail[0]), "tail_fit_const": float(np.exp(p_tail[1])),
            "tail_const_expected_3C": float(3 * C_LIM), "tail_exponent_expected": -2 / 3,
            "k_C_check": {str(int(a)): float(b * a ** (5 / 3) / 1.0) for a, b in zip(kk, cc)},
            "N_asymptotic_tol": {f"{t:g}": float((3 * C_LIM / t) ** 1.5) for t in TOLS}}


# ---- NAFF line table --------------------------------------------------------------------
def _canon(ep):
    return ep.omega * ep.normal


def naff_table(K=24, years=200, dt=10.0, top=10):
    t = MARS.epoch + np.arange(0.0, years * 365.25, dt)
    X = geocentric(MARS, t)
    chain = naff_vector(X, t, K)
    ref = [c for c in geocentric_chain(MARS, EARTH, 40, 40, t_ref=t[0]) if c.omega != 0]
    ref.sort(key=lambda c: -c.radius)
    rows = []
    for e in ref[:top]:
        Om = _canon(e)
        c = min(chain, key=lambda q: np.linalg.norm(_canon(q) - Om))
        Oc = _canon(c)
        n = Om / np.linalg.norm(Om)
        p_ref, p_c = e.position(t[0]), c.position(t[0])
        ph = float(np.arctan2(n @ np.cross(p_ref, p_c), p_ref @ p_c))
        rows.append({"label": f"{e.k:+d} n_{'b' if e.source == 'body' else 'E'}",
                     "source": e.source, "k": e.k,
                     "omega_true": float(e.omega), "omega_recovered": float(c.omega),
                     "period_true_days": float(2 * np.pi / abs(e.omega)),
                     "freq_err_abs": float(np.linalg.norm(Oc - Om)),
                     "freq_err_rel": float(np.linalg.norm(Oc - Om) / np.linalg.norm(Om)),
                     "radius_true_AU": float(e.radius), "radius_recovered_AU": float(c.radius),
                     "radius_err_rel": float(abs(c.radius / e.radius - 1)),
                     "phase_err_rad": ph})
    return {"setup": f"Keplerian geocentric Mars, {years} yr, {dt:g}-day sampling, Hann^4 window, K={K}",
            "reference": "analytic geocentric_chain(N_body=40, N_earth=40), ten largest non-fixed circles",
            "lines": rows}


def print_latex(a, b, naff):
    print("\n% ---- near-parabolic N table (tail = sup_M |z - z_N|/a) ----")
    print(r"$e$ & $N_{\rm cert}(10^{-2})$ & $N_{\rm actual}(10^{-2})$ & $N_{\rm cert}(10^{-6})$ & $N_{\rm actual}(10^{-6})$ \\")
    for r in a:
        cells = []
        for tol in ("0.01", "1e-06"):
            nc = r["N_cert"][tol]
            na = r["N_actual"][tol]
            cells += ["--" if nc is None else f"{nc:,}", f"{na:,}" + ("$^*$" if r["extrapolated"][tol] else "")]
        print(f"{r['e']:g} & " + " & ".join(cells) + r" \\")
    print("% * extrapolated with T ~ 3C N^(-2/3)")
    print(f"% C = {b['C_closed_form']:.5f}; fitted tail exponent {b['tail_fit_exponent']:.4f}, const {b['tail_fit_const']:.4f}")
    print("\n% ---- NAFF line table ----")
    print(r"line & period [d] & $\delta\omega/\omega$ & radius [AU] & $\delta r/r$ & $\delta\phi$ [rad] \\")
    for l in naff["lines"]:
        print(f"${l['label']}$ & {l['period_true_days']:.2f} & {l['freq_err_rel']:.1e} & "
              f"{l['radius_true_AU']:.5f} & {l['radius_err_rel']:.1e} & {l['phase_err_rad']:+.1e} \\\\")


def main():
    a = part_a()
    b = part_b()
    naff = naff_table()
    OUT.mkdir(exist_ok=True)
    (OUT / "near_parabolic.json").write_text(json.dumps(
        {"tolerances": TOLS, "kmax": KMAX, "rows": a, "tail_fit": b}, indent=1))
    (OUT / "naff_lines_mars.json").write_text(json.dumps(naff, indent=1))
    print_latex(a, b, naff)


if __name__ == "__main__":
    main()
