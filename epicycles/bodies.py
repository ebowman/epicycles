"""J2000-epoch element sets for a few bodies.

Planets: JPL SSD "Approximate Positions of the Planets", Table 1 (valid
1800-2050 AD), values at J2000 (JD 2451545.0), rates ignored; fetched
2026-10-14 from https://ssd.jpl.nasa.gov/planets/approx_pos.html. Conversion:
M0 = L - varpi, omega = varpi - Omega. Earth is the Earth-Moon barycenter.
Asteroids: JPL SBDB osculating elements (epoch JD 2461200.5 TDB, heliocentric
ecliptic J2000), fetched 2026-10-14 via https://ssd-api.jpl.nasa.gov/sbdb.api
with full-prec=1. All values fetched, not from memory.
"""
from __future__ import annotations

from math import radians

from .kepler import MU_SUN, Elements

J2000 = 2451545.0


def _planet(name, a, e, I, L, varpi, node, epoch=J2000) -> Elements:
    """Build Elements from a, e [AU, -], I, L, long.peri, long.node [deg]."""
    return Elements(a, e, radians(I), radians(node), radians(varpi - node),
                    radians((L - varpi) % 360.0), epoch, MU_SUN, name)


def _asteroid(name, a, e, I, node, w, ma, epoch=2461200.5) -> Elements:
    """Build Elements from SBDB a, e, i, om, w, ma (degrees; a in AU)."""
    return Elements(a, e, radians(I), radians(node), radians(w), radians(ma),
                    epoch, MU_SUN, name)


EARTH = _planet("Earth", 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0)
MARS = _planet("Mars", 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891)
JUPITER = _planet("Jupiter", 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909)
SATURN = _planet("Saturn", 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448)

CERES = _asteroid("Ceres", 2.765552595034094, 0.07969229514816586, 10.58802780183462,
                  80.24862682043221, 73.29421453021587, 274.4193463761342)
EROS = _asteroid("Eros", 1.458243716760167, 0.2228779627700761, 10.82854410314273,
                 304.2679713350896, 178.9181319135911, 62.51145501986792)
PHAETHON = _asteroid("Phaethon", 1.271464620755913, 0.8896722844020535, 22.31052728999002,
                     265.0988060323288, 322.3001684957205, 301.4858238742286)

BODIES = {b.name.upper(): b for b in
          (EARTH, MARS, JUPITER, SATURN, CERES, EROS, PHAETHON)}
