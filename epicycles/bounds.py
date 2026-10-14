"""Rigorous truncation bounds (derivation Part 4)."""
from __future__ import annotations

import numpy as np


def rho(e: float) -> float:
    """Kapteyn decay ratio e exp(s)/(1+s), s = sqrt(1-e^2)  (derivation 1.7)."""
    s = np.sqrt(1.0 - e * e)
    return float(e * np.exp(s) / (1.0 + s))


def tail_bound(e: float, N: int, a: float = 1.0) -> float:
    """eps_N(a,e) = a * 2(1/e + rho) rho^(N+1) / ((N+1)(1-rho))   (derivation 4.1).

    Bounds sup_t |z - z_N| for the |k|<=N truncation.  e = 0: 0 for N >= 1,
    and a for N = 0 (the single circle is dropped).  Raises ValueError for
    e outside [0,1) or N < 0.
    """
    if not (0.0 <= e < 1.0):
        raise ValueError("eccentricity must satisfy 0 <= e < 1")
    if N < 0:
        raise ValueError("N must be >= 0")
    if e == 0.0:
        return 0.0 if N >= 1 else float(a)
    r = rho(e)
    return float(a * 2.0 * (1.0 / e + r) * r ** (N + 1) / ((N + 1) * (1.0 - r)))


def geocentric_tail_bound(body, earth, N_body: int, N_earth: int) -> float:
    """Derivation (4.2): eps_Nb(body) + eps_Ne(earth) [AU]."""
    return tail_bound(body.e, N_body, body.a) + tail_bound(earth.e, N_earth, earth.a)
