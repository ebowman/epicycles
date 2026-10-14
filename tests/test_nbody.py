import numpy as np

from epicycles.bodies import CERES, EARTH, JUPITER, SATURN
from epicycles.kepler import position
from epicycles.nbody import MASS_RATIOS, energy, integrate, integrate_state

BODIES = {"Earth": EARTH, "Jupiter": JUPITER, "Saturn": SATURN, "Ceres": CERES}
YR = 365.25


def test_energy_conserved():
    # Ceres given a small real-ish mass so that it is part of the conserved energy
    masses = dict(MASS_RATIOS, Ceres=1e-9)
    t, names, m, p, v = integrate_state(BODIES, CERES.epoch, 50 * YR, 10.0, masses)
    E = energy(p, v, m)
    assert np.max(np.abs(E / E[0] - 1)) < 1e-9
    # and the massive subsystem with a massless Ceres
    t, names, m, p, v = integrate_state(BODIES, CERES.epoch, 50 * YR, 10.0)
    E = energy(p, v, m)
    assert np.max(np.abs(E / E[0] - 1)) < 1e-9


def test_massless_planets_match_kepler():
    tr = integrate(BODIES, CERES.epoch, 25 * YR, 5.0, planets_massless=True)
    err = np.max(np.linalg.norm(tr["Ceres"] - position(CERES, tr["t"]), axis=1))
    assert err < 1e-8
    err_e = np.max(np.linalg.norm(tr["Earth"] - position(EARTH, tr["t"]), axis=1))
    assert err_e < 1e-8


def test_perturbed_ceres_deviation_plausible():
    tr = integrate(BODIES, CERES.epoch, 30 * YR, 5.0)
    dev = np.linalg.norm(tr["Ceres"] - position(CERES, tr["t"]), axis=1)
    print("Ceres deviation from Kepler after 30 yr [AU]:", dev[-1], "max", dev.max())
    # osculating vs mean elements: ~0.13 AU (mostly mean-motion offset) after 30 yr
    assert 1e-3 < dev.max() < 1.0
