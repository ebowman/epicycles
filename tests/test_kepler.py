import numpy as np
import pytest

from epicycles.bodies import BODIES, EARTH, MARS, PHAETHON
from epicycles.frames import ecliptic_lonlat, geocentric
from epicycles.kepler import Elements, MU_SUN, mean_motion, position, solve_kepler, velocity


def test_conservation():
    for el in (PHAETHON, MARS, EARTH):
        P = 2 * np.pi / mean_motion(el)
        t = el.epoch + np.linspace(0, 3 * P, 400)
        r, v = position(el, t), velocity(el, t)
        en = 0.5 * np.sum(v**2, -1) - el.mu / np.linalg.norm(r, axis=-1)
        h = np.cross(r, v)
        assert np.ptp(en) / abs(en[0]) < 1e-12
        assert np.max(np.linalg.norm(h - h[0], axis=-1)) / np.linalg.norm(h[0]) < 1e-12


def test_solver_residual():
    M = np.linspace(-4 * np.pi, 4 * np.pi, 2001)
    for e in (0, 0.1, 0.5, 0.9, 0.99):
        E = solve_kepler(M, e)
        assert np.max(np.abs(E - e * np.sin(E) - M)) < 1e-12


def test_circular():
    el = Elements(1.3, 0.0, 0.3, 1.0, 2.0, 0.5, 2451545.0)
    t = el.epoch + np.linspace(0, 500, 50)
    assert np.allclose(np.linalg.norm(position(el, t), axis=-1), 1.3, rtol=1e-14)
    sp = np.linalg.norm(velocity(el, t), axis=-1)
    assert np.allclose(sp, 1.3 * mean_motion(el), rtol=1e-13)


def test_retrograde_mars():
    t0 = 2452880.0
    def rate(t):
        lon, _ = ecliptic_lonlat(geocentric(MARS, np.array([t - 0.5, t + 0.5])))
        d = lon[1] - lon[0]
        return (d + np.pi) % (2 * np.pi) - np.pi
    assert rate(t0) < 0
    assert rate(t0 - 90) > 0
    assert rate(t0 + 90) > 0


@pytest.mark.parametrize("e", [1.0, 1.5, -0.1])
def test_bad_e(e):
    with pytest.raises(ValueError):
        solve_kepler(1.0, e)
    with pytest.raises(ValueError):
        position(Elements(1.0, e, 0, 0, 0, 0, 0.0), 1.0)


def test_wrapping():
    M = np.linspace(-3, 3, 13)
    for e in (0.2, 0.9):
        E1, E2 = solve_kepler(M, e), solve_kepler(M + 2 * np.pi, e)
        d = (E2 - E1 - 2 * np.pi)
        assert np.max(np.abs(d)) < 1e-12


def test_shapes():
    t = MARS.epoch + np.arange(12.0).reshape(3, 4)
    assert position(MARS, t).shape == (3, 4, 3)
    assert velocity(MARS, t).shape == (3, 4, 3)
    assert position(MARS, 2451545.0).shape == (3,)
    assert geocentric(MARS, t).shape == (3, 4, 3)
    lon, lat = ecliptic_lonlat(geocentric(MARS, t))
    assert lon.shape == (3, 4) and np.all((lon >= 0) & (lon < 2 * np.pi))


def test_vis_viva():
    for el in BODIES.values():
        t = el.epoch + np.linspace(0, 2000, 100)
        r = np.linalg.norm(position(el, t), axis=-1)
        v2 = np.sum(velocity(el, t) ** 2, -1)
        assert np.allclose(v2, MU_SUN * (2 / r - 1 / el.a), rtol=1e-11)
