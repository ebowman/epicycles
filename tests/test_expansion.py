"""Tests for the analytic epicycle expansion; names cite derivation.md Part 6 claims."""
import json
from dataclasses import replace

import numpy as np
import pytest
from scipy.special import jv, jvp

from epicycles.bodies import EARTH, MARS, PHAETHON
from epicycles.bounds import geocentric_tail_bound, rho, tail_bound
from epicycles.expansion import (Epicycle, chain_arms, ellipse_to_circles, evaluate_chain,
                                 geocentric_chain, kepler_fourier, to_epicycles)
from epicycles.export import chain_from_json, chain_to_json, write_chain
from epicycles.frames import geocentric
from epicycles.kepler import Elements, MU_SUN, mean_motion, position, solve_kepler

YEAR = 365.25


def _series(e, N, M):
    c = kepler_fourier(e, N)
    k = np.arange(-N, N + 1)
    return np.exp(1j * np.outer(M, k)) @ c


def _true(e, M):
    E = solve_kepler(M, e)
    return (np.cos(E) - e) + 1j * np.sqrt(1 - e * e) * np.sin(E)


def test_claim1_coefficients_vs_direct_solve():
    M = np.linspace(0, 2 * np.pi, 20001)
    for e in (0.0167, 0.1, 0.5):
        assert np.max(np.abs(_series(e, 60, M) - _true(e, M))) < 1e-12


def test_claim2_derivative_forms_agree():
    for e in (0.0167, 0.0934, 0.25, 0.6, 0.89):
        N = 100
        c = kepler_fourier(e, N)
        k = np.arange(1, N + 1)
        s = np.sqrt(1 - e * e)
        ck = (jvp(k, k * e) + s / e * jv(k, k * e)) / k
        cmk = (jvp(k, k * e) - s / e * jv(k, k * e)) / k
        assert np.max(np.abs(c[N + 1:] - ck)) < 1e-14
        assert np.max(np.abs(c[:N][::-1] - cmk)) < 1e-14


def test_claim3_periapsis_sum_rule():
    e = 0.3
    k = np.arange(1, 201)
    assert abs(-1.5 * e + np.sum(2 / k * jvp(k, k * e)) - (1 - e)) < 1e-12


def test_claim4_e0_single_circle():
    c = kepler_fourier(0.0, 5)
    assert c[5 + 1] == 1 and np.count_nonzero(c) == 1
    el = Elements(1.3, 0.0, 0.3, 1.0, 2.0, 0.5, 2451545.0)
    ch = to_epicycles(el, 10)
    assert len(ch) == 1 and ch[0].radius == 1.3 and ch[0].k == 1
    t = el.epoch + np.linspace(0, 400, 50)
    assert np.allclose(evaluate_chain(ch, t), position(el, t), atol=1e-13)


def test_e0_continuity_and_tiny_e_stability():
    c0 = kepler_fourier(0.0, 6)
    c1 = kepler_fourier(1e-10, 6)
    assert np.max(np.abs(c1 - c0)) < 1e-9
    assert abs(c1[6 + 1] - 1) < 1e-14
    assert abs(c1[6 - 1] / 1e-20 - 1 / 8) < 1e-6  # c_{-1} ~ e^2/8, no cancellation damage
    assert abs(c1[6 + 2] / 1e-10 - 0.5) < 1e-6


def test_kepler_fourier_errors():
    for bad in (-0.1, 1.0, 1.5):
        with pytest.raises(ValueError):
            kepler_fourier(bad, 3)
    with pytest.raises(ValueError):
        kepler_fourier(0.5, -1)


def test_claim5_small_e_coefficients():
    e = 1e-3
    c = kepler_fourier(e, 3)
    assert c[3] == -1.5 * e
    assert abs(c[3 + 2] / e - 0.5) < 1e-3 * 0.5
    assert abs(c[3 - 1] / e**2 - 0.125) < 1e-3 * 0.125
    assert abs((1 - c[3 + 1]) / (e**2 / 2) - 1) < 1e-3


def test_claim6_kapteyn_and_coefficient_bounds():
    k = np.arange(1, 400)
    for e in (0.0167, 0.0934, 0.25, 0.6, 0.89, 0.99):
        r = rho(e)
        assert np.all(np.abs(jv(k, k * e)) <= r**k + 1e-15)
        c = kepler_fourier(e, 399)
        s = np.sqrt(1 - e * e)
        assert np.all(np.abs(c[399 + 1:]) <= r**k / k * ((1 + s) / e + r) + 1e-15)
        assert np.all(np.abs(c[:399][::-1]) <= r**k / k * ((1 - s) / e + r) + 1e-15)


@pytest.mark.parametrize("e", [0.0167, 0.0934, 0.6, 0.89])
def test_claim7_truncation_bound_holds(e):
    M = np.linspace(0, np.pi, 3001)  # error symmetric under M -> -M
    Ns = [0, 1, 2, 3, 5, 8, 13, 21, 40, 80, 150, 300] + ([600, 1000] if e > 0.8 else [])
    zt = _true(e, M)
    for N in Ns:
        err = np.max(np.abs(_series(e, N, M) - zt))
        assert err <= tail_bound(e, N, 1.0) + 1e-13, (e, N)


def test_tail_bound_edge_cases():
    assert tail_bound(0.0, 3, 2.0) == 0.0 and tail_bound(0.0, 0, 2.0) == 2.0
    with pytest.raises(ValueError):
        tail_bound(1.0, 3)
    with pytest.raises(ValueError):
        tail_bound(0.3, -1)


def _rand_el(rng, e=0.5):
    return Elements(1.7, e, rng.uniform(0, 3), rng.uniform(0, 6.28), rng.uniform(0, 6.28),
                    rng.uniform(0, 6.28), 2451545.0)


def test_claim9_3d_rotation_and_to_epicycles():
    rng = np.random.default_rng(1)
    for _ in range(3):
        el = _rand_el(rng)
        ch = to_epicycles(el, 60)
        for ep in ch:
            assert np.allclose(np.cross(ep.u, ep.v), ep.normal, atol=1e-14)
            assert abs(ep.u @ ep.v) < 1e-14 and abs(np.linalg.norm(ep.normal) - 1) < 1e-14
        t = el.epoch + np.linspace(0, 3 * 2 * np.pi / mean_motion(el), 2000)
        assert np.max(np.abs(evaluate_chain(ch, t) - position(el, t))) < 1e-12 * el.a
        # t_ref independence
        ch2 = to_epicycles(el, 60, t_ref=el.epoch + 123.4)
        assert np.max(np.abs(evaluate_chain(ch2, t) - evaluate_chain(ch, t))) < 1e-12


def test_claim10_two_circle_decomposition():
    rng = np.random.default_rng(2)
    for _ in range(5):
        el = _rand_el(rng)
        k = int(rng.integers(1, 11))
        c = kepler_fourier(el.e, k)
        ch = to_epicycles(el, k)
        P, Q, Nn = ch[0].u, ch[0].v, ch[0].normal
        A = el.a * (c[k + k] * (P - 1j * Q) + c[k - k] * (P + 1j * Q))
        circles = ellipse_to_circles(A, k * mean_motion(el), t_ref=0.0)
        M = np.linspace(0, 6, 40)
        t = M / mean_motion(el)  # small t: avoids JD rounding
        want = np.real(A[None] * np.exp(1j * k * M)[:, None])
        assert np.max(np.abs(evaluate_chain(circles, t) - want)) < 1e-14
        for ep in circles:  # every circle lies in the orbital plane (normal = +-N)
            assert abs(ep.normal @ Nn) > 1 - 1e-12 and abs(ep.normal @ P) < 1e-12
        assert len(circles) == 2
        rr = sorted(ep.radius for ep in circles)
        ck, cmk = abs(c[2 * k]), abs(c[0])
        assert abs(rr[1] + rr[0] - el.a * (ck + cmk)) < 1e-14  # semi-major
        sv = np.linalg.svd(np.stack([A.real, -A.imag], 1), compute_uv=False)
        assert abs(sv[0] - el.a * (ck + cmk)) < 1e-13 and abs(sv[1] - el.a * abs(ck - cmk)) < 1e-13


def test_ellipse_to_circles_random_A():
    rng = np.random.default_rng(3)
    t = np.linspace(-5, 5, 200)
    for _ in range(20):
        A = rng.normal(size=3) + 1j * rng.normal(size=3)
        w = rng.uniform(-2, 2)
        ch = ellipse_to_circles(A, w, t_ref=0.3)
        want = np.real(A[None] * np.exp(1j * w * (t - 0.3))[:, None])
        assert np.max(np.abs(evaluate_chain(ch, t) - want)) < 1e-13
        # semi-axes = singular values of [Re A, -Im A]
        sv = np.linalg.svd(np.stack([A.real, -A.imag], 1), compute_uv=False)
        r = [c.radius for c in ch]
        assert abs((r[0] + r[1]) - sv[0]) < 1e-12 and abs(abs(r[0] - r[1]) - sv[1]) < 1e-12
    # degenerate (linear) and zero cases
    lin = ellipse_to_circles(np.array([1.0, 2.0, 0.5]) * (0.3 + 0.4j), 1.1)
    want = np.real((np.array([1.0, 2.0, 0.5]) * (0.3 + 0.4j))[None] * np.exp(1j * 1.1 * t)[:, None])
    assert np.max(np.abs(evaluate_chain(lin, t) - want)) < 1e-13
    assert ellipse_to_circles(np.zeros(3, complex), 1.0) == []


def test_chain_sorted_arms_and_negation():
    ch = geocentric_chain(MARS, EARTH, 10, 10)
    r = [c.radius for c in ch]
    assert r == sorted(r, reverse=True)
    assert {c.source for c in ch} == {"body", "earth"}
    t = MARS.epoch + np.array([0.0, 100.0, 5000.0])
    arms = chain_arms(ch, t)
    assert arms.shape == (3, len(ch) + 1, 3) and np.all(arms[:, 0] == 0)
    assert np.allclose(arms[:, -1], evaluate_chain(ch, t), atol=1e-13)
    assert chain_arms(ch, MARS.epoch).shape == (len(ch) + 1, 3)
    assert evaluate_chain(ch, MARS.epoch).shape == (3,)


def _geo_err(body, Nb, Ne, years, n=30001):
    t = body.epoch + np.linspace(0, years * YEAR, n)
    return np.max(np.abs(evaluate_chain(geocentric_chain(body, EARTH, Nb, Ne), t) - geocentric(body, t)))


def test_claim11_mars_convergence_and_bound():
    errs = [_geo_err(MARS, N, N, 60) for N in range(1, 13)]
    # measured: monotone decreasing from N=1 (derivation notes non-monotonicity only
    # for planar |k|<=2 subsets in longitude); floor ~1e-13 is reached at N~13.
    assert all(b < a for a, b in zip(errs, errs[1:])), errs
    for N, err in zip(range(1, 11), errs):
        assert err <= geocentric_tail_bound(MARS, EARTH, N, N) + 1e-13


def test_mars_n40_accuracy():
    assert _geo_err(MARS, 40, 40, 60) < 1e-10


def test_phaethon_convergence():
    # Task target was <1e-6 AU at N=200; MEASURED 4.24e-6 (FLAGGED). 1e-6 needs N~235.
    e200 = _geo_err(PHAETHON, 200, 40, 60, 40001)
    e240 = _geo_err(PHAETHON, 240, 40, 60, 40001)
    assert e200 < 5e-6
    assert e240 < 1e-6
    assert e200 <= geocentric_tail_bound(PHAETHON, EARTH, 200, 40)


def _spectrum_lines(b_el, years=200.0, dt=2.0):
    t = b_el.epoch + np.arange(0, years * YEAR, dt)
    r = geocentric(b_el, t)
    # sin^8 window (Laskar's Hann^p, p=4): far-sidelobe decay ~m^-9.  Plain Hann
    # (p=1) leaks 1e-4 at 12 bins and cannot resolve 1e-6 near strong lines.
    w = np.sin(np.pi * (np.arange(len(t)) + 0.5) / len(t)) ** 8
    F = np.sqrt(sum(np.abs(np.fft.rfft((r[:, i] - r[:, i].mean()) * w)) ** 2 for i in range(3)))
    f = np.fft.rfftfreq(len(t), dt)  # cycles/day
    return f, F / F.max()


def test_claim13_no_cross_terms_in_cartesian_spectrum():
    """Measured (see summary): out-of-mask relative amplitude vs 1e-6 proposal."""
    f, F = _spectrum_lines(MARS)
    fb, fe = mean_motion(MARS) / (2 * np.pi), mean_motion(EARTH) / (2 * np.pi)
    df = f[1]
    lines = np.array([k * fb for k in range(1, 41)] + [j * fe for j in range(1, 41)])
    dist = np.min(np.abs(f[:, None] - lines[None]), axis=1)
    out = F[(dist > 12 * df) & (f > 0.2 * fb)]
    assert out.max() < 1e-6, out.max()
    # explicit combination frequencies |k fb +- j fe|, away from any base line
    checked = 0
    for k in (1, 2, 3):
        for j in (1, 2):
            for sgn in (1, -1):
                fc = abs(k * fb + sgn * j * fe)
                if np.min(np.abs(lines - fc)) > 13 * df:
                    i = np.argmin(np.abs(f - fc))
                    assert F[i - 2:i + 3].max() < 1e-6
                    checked += 1
    assert checked >= 3
    # sanity: the base lines are really visible
    i = np.argmin(np.abs(f - fe))
    assert F[i - 3:i + 4].max() > 1e-2


def _synthetic(ratio):
    """Body with n_b/n_E = ratio on an inclined orbit, same mu."""
    a_b = EARTH.a * ratio ** (-2 / 3)
    return Elements(a_b, 0.2, 0.2, 0.5, 1.0, 0.3, 0.0, MU_SUN, "synthetic")


def test_claim14_rational_ratio_periodic_irrational_not():
    """Rational n_b/n_E = 1/2: exact periodicity; 1/sqrt(2): no near-recurrence.

    Epochs are 0 (not J2000 JD) so t has no 1e-10 d rounding; see summary for values.
    """
    E0 = replace(EARTH, epoch=0.0)
    Tp = 4 * np.pi / mean_motion(E0)
    body = _synthetic(0.5)
    assert abs(mean_motion(body) / mean_motion(E0) - 0.5) < 1e-12
    t = np.linspace(0, 6 * Tp, 3000)
    assert np.max(np.abs(geocentric(body, t + Tp, E0) - geocentric(body, t, E0))) < 1e-12
    irr = _synthetic(1 / np.sqrt(2))
    tw = np.linspace(0, 12 * Tp, 600)
    base = geocentric(irr, tw, E0)
    tps = np.linspace(Tp, 10 * Tp, 4000)
    m = min(np.max(np.abs(geocentric(irr, tw + tp, E0) - base)) for tp in tps)
    assert m > 1e-3  # bounded away from 0 (amplitude scale ~ 1-2 AU)
    assert np.max(np.abs(geocentric(irr, tw + Tp, E0) - base)) > 1e-1


# ---- export -------------------------------------------------------------------

def test_export_round_trip(tmp_path):
    ch = geocentric_chain(MARS, EARTH, 12, 6)
    meta = {"body": "Mars", "earth": "Earth", "N_body": 12, "N_earth": 6, "t_ref": MARS.epoch + 50.0}
    p = tmp_path / "sub" / "mars.json"
    write_chain(p, ch, meta)
    d = json.loads(p.read_text())
    assert d["schema_version"] == 1 and d["meta"]["units"] == {"length": "AU", "time": "day", "angle": "rad"}
    assert d["meta"]["t_ref"] == MARS.epoch + 50.0 and len(d["epicycles"]) == len(ch)
    assert set(d["epicycles"][0]) == {"radius", "omega", "phase", "normal", "u", "v", "source", "k"}
    t = MARS.epoch + np.linspace(0, 1000, 50)
    # evaluate straight from the documented formula, independent of the Python classes
    pos = np.zeros((50, 3))
    for e in d["epicycles"]:
        th = e["phase"] + e["omega"] * (t - d["meta"]["t_ref"])
        pos += e["radius"] * (np.cos(th)[:, None] * np.array(e["u"]) + np.sin(th)[:, None] * np.array(e["v"]))
    assert np.max(np.abs(pos - evaluate_chain(ch, t))) < 1e-12
    assert np.max(np.abs(evaluate_chain(chain_from_json(d), t) - evaluate_chain(ch, t))) < 1e-12
    assert chain_to_json(ch, meta) == d
