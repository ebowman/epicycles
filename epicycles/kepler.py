"""Kepler propagation in the ecliptic J2000 frame.

Conventions (match theory/derivation.md): elements a, e, i, Omega, omega, M0,
angles in radians; M = M0 + n (t - epoch); n = sqrt(mu / a^3); perifocal x axis
toward periapsis; R = R_z(Omega) R_x(i) R_z(omega). Units: AU, days.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Gaussian gravitational constant, rad/day (sqrt of mu_sun in AU^3/day^2 units).
GAUSS_K = 0.01720209895
#: Solar gravitational parameter, AU^3 / day^2.
MU_SUN = GAUSS_K**2


@dataclass(frozen=True)
class Elements:
    """Osculating Keplerian elements.

    a: semi-major axis [AU]; e: eccentricity (0 <= e < 1); i, Omega, omega:
    inclination, node longitude, argument of periapsis [rad]; M0: mean anomaly
    at ``epoch`` [rad]; epoch: Julian date (TDB) [days]; mu: gravitational
    parameter [AU^3/day^2]; name: label.
    """

    a: float
    e: float
    i: float
    Omega: float
    omega: float
    M0: float
    epoch: float
    mu: float = MU_SUN
    name: str = ""


def solve_kepler(M, e, tol: float = 1e-15, maxiter: int = 100):
    """Solve Kepler's equation M = E - e sin E for E.

    M: mean anomaly [rad], scalar or array (any value; wrapping handled, and
    E - M is periodic in M). e: eccentricity, scalar or broadcastable, 0 <= e < 1.
    Returns E [rad] (same broadcast shape; float for scalar input), with
    |E - e sin E - M| ~ 1e-15 .. 1e-14. Raises ValueError if e outside [0, 1).
    """
    M = np.asarray(M, dtype=float)
    e = np.asarray(e, dtype=float)
    if np.any(~(e >= 0.0)) or np.any(e >= 1.0):
        raise ValueError("eccentricity must satisfy 0 <= e < 1")
    M, e = np.broadcast_arrays(M, e)
    Mr = np.remainder(M + np.pi, 2 * np.pi) - np.pi  # in [-pi, pi)
    E = np.where(e > 0.8, np.pi * np.sign(Mr) + 0 * Mr, Mr + e * np.sin(Mr))
    E = np.where((e > 0.8) & (Mr == 0), 0.0, E)
    for _ in range(maxiter):
        f = E - e * np.sin(E) - Mr
        dE = f / (1.0 - e * np.cos(E))
        dE = np.clip(dE, -1.0, 1.0)
        E = E - dE
        if np.all(np.abs(dE) < tol) or np.all(np.abs(E - e * np.sin(E) - Mr) < 1e-15):
            break
    E = E + (M - Mr)
    return float(E) if E.ndim == 0 else E


def mean_motion(el: Elements) -> float:
    """Mean motion n = sqrt(mu / a^3) [rad/day]."""
    return float(np.sqrt(el.mu / el.a**3))


def _rotation(el: Elements) -> np.ndarray:
    """Rotation matrix R_z(Omega) R_x(i) R_z(omega): perifocal -> ecliptic (3x3)."""
    def rz(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

    c, s = np.cos(el.i), np.sin(el.i)
    rx = np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])
    return rz(el.Omega) @ rx @ rz(el.omega)


def _check(el: Elements) -> None:
    if not (0.0 <= el.e < 1.0):
        raise ValueError("eccentricity must satisfy 0 <= e < 1")
    if not el.a > 0:
        raise ValueError("semi-major axis must be positive")


def _anomaly(el: Elements, t):
    t = np.asarray(t, dtype=float)
    M = el.M0 + mean_motion(el) * (t - el.epoch)
    return solve_kepler(M, el.e)


def position(el: Elements, t) -> np.ndarray:
    """Heliocentric ecliptic J2000 position [AU].

    t: Julian date TDB [days], scalar or array of shape (...). Returns (...,3).
    """
    _check(el)
    E = np.asarray(_anomaly(el, t))
    s = np.sqrt(1.0 - el.e**2)
    pf = np.stack([el.a * (np.cos(E) - el.e), el.a * s * np.sin(E), np.zeros_like(E)], axis=-1)
    return pf @ _rotation(el).T


def velocity(el: Elements, t) -> np.ndarray:
    """Heliocentric ecliptic J2000 velocity [AU/day].

    t: Julian date TDB [days], scalar or array of shape (...). Returns (...,3).
    """
    _check(el)
    E = np.asarray(_anomaly(el, t))
    n = mean_motion(el)
    s = np.sqrt(1.0 - el.e**2)
    f = n * el.a / (1.0 - el.e * np.cos(E))  # dE/dt * a
    vf = np.stack([-f * np.sin(E), f * s * np.cos(E), np.zeros_like(E)], axis=-1)
    return vf @ _rotation(el).T
