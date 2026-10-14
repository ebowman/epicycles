"""Analytic epicycle expansion of Kepler orbits (theory/derivation.md, Parts 1-3).

Conventions
-----------
* Units: AU, days, radians.
* An :class:`Epicycle` is one uniformly rotating circle (a "term") in a plane
  with unit normal ``normal`` and orthonormal in-plane basis ``u, v`` with
  ``u x v = normal``.  Its position contribution at Julian date ``t`` is::

      radius * (cos(theta) * u + sin(theta) * v),
      theta = phase + omega * (t - t_ref)

  ``omega`` is signed: omega > 0 rotates counter-clockwise about +normal,
  omega < 0 clockwise.  ``phase`` is the angle from ``u`` at ``t = t_ref``.
* The full motion is the *sum* of all terms (all circles are concentric on the
  origin; the nested-arm picture is the same sum drawn tip to tail, see
  :func:`chain_arms`).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.special import jv

from .kepler import Elements, _check, _rotation, mean_motion


@dataclass(frozen=True, eq=False)
class Epicycle:
    """One circle of an epicycle chain; see module docstring for conventions.

    radius [AU] >= 0; omega [rad/day] signed; phase [rad] at ``t_ref`` [JD];
    normal, u, v: unit 3-vectors (ecliptic J2000), u x v = normal.
    source: "body" / "earth" / "" ; k: harmonic index (None if not applicable).
    """

    radius: float
    omega: float
    phase: float
    normal: np.ndarray
    u: np.ndarray
    v: np.ndarray
    t_ref: float = 0.0
    source: str = ""
    k: int | None = None

    def position(self, t) -> np.ndarray:
        """Contribution of this circle at time(s) t, shape (...,3)."""
        t = np.asarray(t, dtype=float)
        th = self.phase + self.omega * (t - self.t_ref)
        return (self.radius * np.cos(th))[..., None] * self.u + (
            self.radius * np.sin(th))[..., None] * self.v

    def rereferenced(self, t_ref: float) -> "Epicycle":
        """Same circle with phase expressed at a different reference epoch."""
        return replace(self, phase=self.phase + self.omega * (t_ref - self.t_ref), t_ref=t_ref)

    def negated(self) -> "Epicycle":
        """The circle traversed with opposite sign (phase + pi)."""
        return replace(self, phase=self.phase + np.pi)


def kepler_fourier(e: float, N: int) -> np.ndarray:
    """Fourier coefficients c_k, k = -N..N, of z/a = sum_k c_k exp(i k M).

    Returns a real-valued complex array of length 2N+1; index ``k + N`` holds
    c_k (so ``c[N]`` is c_0).  z is the perifocal position x + iy (focus at the
    origin, x to periapsis).  Formulas (derivation 1.4, 1.5):
    c_0 = -3e/2; c_k = (1/k)[(1+s)/e J_k(ke) - J_{k+1}(ke)];
    c_-k = (1/k)[(1-s)/e J_k(ke) - J_{k+1}(ke)], s = sqrt(1-e^2).
    For stability (1-s) is evaluated as e^2/(1+s), so (1-s)/e = e/(1+s), and
    J_k(ke)/e is formed from jv directly (relatively accurate for tiny e).
    e = 0 is special-cased exactly (c_1 = 1, all others 0).
    Raises ValueError for e < 0, e >= 1 or N < 0.
    """
    if not (0.0 <= e < 1.0):
        raise ValueError("eccentricity must satisfy 0 <= e < 1")
    if N < 0 or int(N) != N:
        raise ValueError("N must be a non-negative integer")
    N = int(N)
    c = np.zeros(2 * N + 1, dtype=complex)
    if e == 0.0:
        if N >= 1:
            c[N + 1] = 1.0
        return c
    s = np.sqrt(1.0 - e * e)
    c[N] = -1.5 * e
    k = np.arange(1, N + 1)
    jk = jv(k, k * e)
    jk1 = jv(k + 1, k * e)
    c[N + 1:] = ((1.0 + s) / e * jk - jk1) / k
    c[:N][::-1] = (e / (1.0 + s) * jk - jk1) / k
    return c


def _plane_basis(el: Elements):
    R = _rotation(el)
    return R[:, 0].copy(), R[:, 1].copy(), R[:, 2].copy()  # P, Q, N


def to_epicycles(el: Elements, N: int, t_ref: float | None = None,
                 source: str = "body") -> list[Epicycle]:
    """Truncated (|k| <= N) epicycle expansion of the heliocentric orbit.

    Each term a*c_k*exp(i k M) becomes one circle in the orbital plane with
    u = P, v = Q, normal = N (derivation Part 2.1/2.5), radius a|c_k|, signed
    omega = k n, phase = k M(t_ref) + arg c_k (arg c_k in {0, pi}).  The k = 0
    term is a fixed arm (omega = 0, phase = pi, radius 3ae/2).  Terms with
    exactly zero radius are omitted (so e = 0 gives a single circle).
    t_ref defaults to el.epoch.  Sorted by decreasing radius.
    """
    _check(el)
    if t_ref is None:
        t_ref = el.epoch
    c = kepler_fourier(el.e, N)
    P, Q, Nn = _plane_basis(el)
    n = mean_motion(el)
    M_ref = el.M0 + n * (t_ref - el.epoch)
    out = []
    for idx, ck in enumerate(c):
        k = idx - N
        if ck == 0:
            continue
        out.append(Epicycle(el.a * abs(ck.real), k * n,
                            k * M_ref + (np.pi if ck.real < 0 else 0.0),
                            Nn.copy(), P.copy(), Q.copy(), t_ref, source, k))
    out.sort(key=lambda ep: -ep.radius)
    return out


def ellipse_to_circles(A, omega: float, t_ref: float = 0.0, source: str = "",
                       k: int | None = None) -> list[Epicycle]:
    """Decompose Re(A exp(i omega (t - t_ref))), A complex 3-vector, into circles.

    The motion is an ellipse in the plane spanned by Re A and Im A; it equals a
    prograde circle (rate +omega about the plane normal) plus a retrograde one
    (rate -omega) (derivation 2.4).  With plane basis (u, v, n = u x v):
    alpha = (A.u + i A.v)/2, beta = (A.u - i A.v)/2 (plain dot products); then
    A = alpha (u - i v) + beta (u + i v), the prograde circle has radius |alpha|,
    phase arg alpha, rate +omega; the retrograde has radius |beta|, phase
    -arg beta, rate -omega.  Zero-radius circles are omitted.  A degenerate
    (linear) ellipse gets an arbitrary in-plane basis containing its line.
    """
    A = np.asarray(A, dtype=complex)
    Ar, Ai = A.real, A.imag
    scale = max(np.linalg.norm(Ar), np.linalg.norm(Ai))
    if scale == 0.0:
        return []
    nrm = np.cross(Ar, Ai)
    if np.linalg.norm(nrm) > 1e-12 * scale**2:
        n = nrm / np.linalg.norm(nrm)
        u = Ar / np.linalg.norm(Ar) if np.linalg.norm(Ar) > 0 else Ai / np.linalg.norm(Ai)
    else:  # linear motion along w
        w = Ar if np.linalg.norm(Ar) >= np.linalg.norm(Ai) else Ai
        u = w / np.linalg.norm(w)
        helper = np.eye(3)[np.argmin(np.abs(u))]
        n = np.cross(u, helper)
        n /= np.linalg.norm(n)
    v = np.cross(n, u)
    alpha = (A @ u + 1j * (A @ v)) / 2
    beta = (A @ u - 1j * (A @ v)) / 2
    out = []
    if abs(alpha) > 0:
        out.append(Epicycle(abs(alpha), omega, float(np.angle(alpha)), n, u, v, t_ref, source, k))
    if abs(beta) > 0:
        out.append(Epicycle(abs(beta), -omega, -float(np.angle(beta)), n, u, v, t_ref, source, k))
    out.sort(key=lambda ep: -ep.radius)
    return out


def geocentric_chain(body: Elements, earth: Elements, N_body: int, N_earth: int,
                     t_ref: float | None = None) -> list[Epicycle]:
    """Epicycles of r_body - r_earth: body terms plus negated Earth terms.

    Negation is phase + pi (derivation 3.1: no cross terms).  All phases are
    referenced to a common t_ref (default body.epoch).  Sorted by radius desc.
    """
    if t_ref is None:
        t_ref = body.epoch
    chain = to_epicycles(body, N_body, t_ref, "body") + [
        ep.negated() for ep in to_epicycles(earth, N_earth, t_ref, "earth")]
    chain.sort(key=lambda ep: -ep.radius)
    return chain


def _stack(chain):
    return (np.array([c.radius for c in chain]), np.array([c.omega for c in chain]),
            np.array([c.phase for c in chain]), np.array([c.t_ref for c in chain]),
            np.array([c.u for c in chain]).reshape(-1, 3),
            np.array([c.v for c in chain]).reshape(-1, 3))


def evaluate_chain(chain: list[Epicycle], t) -> np.ndarray:
    """Sum of all circles at time(s) t (JD): array shape (...,3) [AU]."""
    t = np.asarray(t, dtype=float)
    if not chain:
        return np.zeros(t.shape + (3,))
    R, w, ph, tr, U, V = _stack(chain)
    th = ph + w * (t[..., None] - tr)
    return (R * np.cos(th)) @ U + (R * np.sin(th)) @ V


def chain_arms(chain: list[Epicycle], t) -> np.ndarray:
    """Cumulative arm tips for drawing, shape (..., len(chain)+1, 3).

    Point 0 is the origin; point j is the sum of the first j circles (in chain
    order); the last point equals evaluate_chain(chain, t).
    """
    t = np.asarray(t, dtype=float)
    parts = np.stack([c.position(t) for c in chain], axis=-2) if chain else np.zeros(t.shape + (0, 3))
    cum = np.cumsum(parts, axis=-2)
    return np.concatenate([np.zeros(t.shape + (1, 3)), cum], axis=-2)
