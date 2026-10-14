"""Geocentric positions and ecliptic coordinates."""
from __future__ import annotations

import numpy as np

from .bodies import EARTH
from .kepler import Elements, position


def geocentric(body: Elements, t, earth: Elements = EARTH) -> np.ndarray:
    """Geocentric ecliptic J2000 position of `body` [AU] = r_body - r_earth.

    t: Julian date TDB [days], scalar or array (...). Returns (...,3).
    """
    return position(body, t) - position(earth, t)


def ecliptic_lonlat(xyz) -> tuple[np.ndarray, np.ndarray]:
    """Ecliptic longitude in [0, 2*pi) and latitude [rad] of vector(s) xyz (...,3)."""
    xyz = np.asarray(xyz, dtype=float)
    lon = np.mod(np.arctan2(xyz[..., 1], xyz[..., 0]), 2 * np.pi)
    lat = np.arctan2(xyz[..., 2], np.hypot(xyz[..., 0], xyz[..., 1]))
    return lon, lat
