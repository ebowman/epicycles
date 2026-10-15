"""Tests for NAFF extraction (epicycles/naff.py; derivation Part 5).

Measured numbers quoted in comments are from the runs that fixed the tolerances.
"""
import time

import numpy as np
import pytest

from epicycles.bodies import EARTH, MARS
from epicycles.expansion import evaluate_chain, geocentric_chain
from epicycles.frames import geocentric
from epicycles.naff import hann_window, naff_complex, naff_vector

W0 = np.array([1.0, -0.7311, 2.3719, 0.4163, -1.6180])          # incommensurate rad/unit
A0 = np.array([1.0, 0.6 * np.exp(.3j), 0.35 * np.exp(-1j), 0.2 * np.exp(2j), 0.08 * np.exp(.5j)])
SLOW = 0.4163
T_SLOW = 2 * np.pi / SLOW


def _synth(n_periods=100, N=4000, t0=5.0):
    T = n_periods * T_SLOW
    t = np.linspace(0, T, N, endpoint=False) + t0
    z = sum(a * np.exp(1j * w * (t - t[0])) for w, a in zip(W0, A0))
    return t, z, T


def _compare(res, w0=W0, a0=A0):
    """Match extracted lines to truth by amplitude order; return (max rel freq err, max rel amp err)."""
    res = sorted(res, key=lambda r: -abs(r[1]))
    order = np.argsort(-abs(a0))
    fe = max(abs(r[0] - w0[i]) / abs(w0[i]) for r, i in zip(res, order))
    ae = max(abs(r[1] - a0[i]) / abs(a0[i]) for r, i in zip(res, order))
    return fe, ae


def test_window_normalised_and_symmetric():
    w = hann_window(101, 4)
    assert w.sum() == pytest.approx(1.0) and np.allclose(w, w[::-1]) and w[0] == 0.0


def test_synthetic_five_exponentials_100_periods():
    t, z, _ = _synth()
    res = naff_complex(z, t, 5)
    assert len(res) == 5
    fe, ae = _compare(res)          # measured: 1.5e-10, 5.7e-8 (p=2..4 all < 2.7e-10, 1.4e-7)
    assert fe < 1e-8 and ae < 1e-6


def test_noise_robustness_40dB():
    # complex white noise, per-sample rms(noise) = rms(signal)/100  (SNR 40 dB)
    t, z, _ = _synth()
    rng = np.random.default_rng(1)
    sig = np.sqrt(np.mean(abs(z) ** 2)) / 100
    zn = z + sig * (rng.normal(size=z.size) + 1j * rng.normal(size=z.size)) / np.sqrt(2)
    fe, ae = _compare(naff_complex(zn, t, 5))   # measured: 1.3e-5, 1.7e-2 (weakest line, amp 0.08)
    assert fe < 1e-4 and ae < 0.05


def test_noise_20dB_still_finds_frequencies():
    t, z, _ = _synth()
    rng = np.random.default_rng(2)
    sig = np.sqrt(np.mean(abs(z) ** 2)) / 10
    zn = z + sig * (rng.normal(size=z.size) + 1j * rng.normal(size=z.size)) / np.sqrt(2)
    fe, ae = _compare(naff_complex(zn, t, 5))   # measured: 4.9e-5, 7.9e-2
    assert fe < 1e-3 and ae < 0.3


def test_real_valued_input():
    t, _, _ = _synth()
    x = np.cos(SLOW * (t - t[0]) + .3) + 0.5 * np.sin(1.7 * (t - t[0]))   # real
    res = naff_complex(x, t, 4)
    om = sorted(r[0] for r in res)
    assert np.allclose(om, [-1.7, -SLOW, SLOW, 1.7], rtol=1e-8)
    assert np.allclose(sorted(abs(r[1]) for r in res), [0.25, 0.25, 0.5, 0.5], rtol=1e-6)


def test_vector_real_input_ellipse_and_constant():
    t = np.linspace(0, 80 * T_SLOW, 3000, endpoint=False) + 2.45e6
    tt = t - t[0]
    xyz = np.stack([1.0 * np.cos(SLOW * tt) + 0.3, 0.4 * np.sin(SLOW * tt + 1) - 0.2,
                    0.1 * np.cos(2.7 * tt) + 0.05], axis=1)
    ch = naff_vector(xyz, t, 2)
    assert np.max(np.abs(evaluate_chain(ch, t) - xyz)) < 1e-6
    fixed = [c for c in ch if c.omega == 0]
    assert len(fixed) == 1 and fixed[0].radius == pytest.approx(np.linalg.norm([0.3, -0.2, 0.05]), rel=1e-8)
    assert all(c.t_ref == t[0] for c in ch)


@pytest.mark.parametrize("sep_bins", [0.1, 0.5, 1.0, 2.0])
def test_near_degenerate_lines_graceful(sep_bins):
    # Two lines closer than ~(p+1) = 5 bins (bin = 2 pi/T) are not resolved with
    # the default Hann^4 window: measured, they come back merged or as a
    # badly-conditioned pair with wrong amplitudes/frequencies.  Required here:
    # no crash and finite output.
    t, _, T = _synth()
    dw = sep_bins * 2 * np.pi / T
    z = np.exp(1j * (t - t[0])) + 0.5 * np.exp(1j * (1.0 + dw) * (t - t[0]))
    res = naff_complex(z, t, 2)
    assert 1 <= len(res) <= 2
    assert all(np.isfinite(w) and np.isfinite(a) for w, a in res)


def test_resolution_limit_8_bins_resolved():
    # measured: separation 4 bins -> freq error 0.012 bin, amp err 0.4%; 8 bins -> < 1e-3 bin
    t, _, T = _synth()
    dw = 8 * 2 * np.pi / T
    z = np.exp(1j * (t - t[0])) + 0.5 * np.exp(1j * (1.0 + dw) * (t - t[0]))
    res = sorted(naff_complex(z, t, 2), key=lambda r: -abs(r[1]))
    assert abs(res[0][0] - 1.0) < 1e-3 * dw / 8 and abs(res[1][0] - (1 + dw)) < 1e-3 * dw / 8
    assert abs(abs(res[1][1]) - 0.5) < 1e-4


@pytest.mark.parametrize("N", [8, 12, 16, 40])
def test_short_windows_finite(N):
    tt = np.arange(N) * 1.0          # 1.5 - 6 periods... even < 1 period of the slow line
    zz = np.exp(1j * 0.3 * tt) + 0.3 * np.exp(-1j * 1.1 * tt)
    res = naff_complex(zz, tt, 3)
    assert all(np.isfinite(w) and np.isfinite(a) for w, a in res)
    ch = naff_vector(np.stack([np.cos(.3 * tt), np.sin(.3 * tt), 0 * tt], 1), tt, 3)
    assert all(np.isfinite(c.radius) and np.isfinite(c.omega) and np.isfinite(c.phase) for c in ch)
    if N == 40:   # measured: ~1.9 periods of the 0.3 line, still 1.7e-7 relative frequency
        assert abs(res[0][0] - 0.3) < 1e-5


def test_input_validation():
    with pytest.raises(ValueError):
        naff_complex(np.ones(5, complex), np.arange(5.0), 1)
    t = np.arange(20.0); t[7] += 0.3
    with pytest.raises(ValueError):
        naff_complex(np.ones(20, complex), t, 1)
    with pytest.raises(ValueError):
        naff_vector(np.ones((20, 2)), np.arange(20.0), 1)


# ---- Kepler geocentric trajectory --------------------------------------------------

def _canon(ep):
    """Angular-velocity vector omega*normal: a circle is (omega, n) ~ (-omega, -n)."""
    return ep.omega * ep.normal


def test_mars_geocentric_matches_analytic_chain():
    t = MARS.epoch + np.arange(0.0, 200 * 365.25, 10.0)        # 200 yr, 10-day sampling
    X = geocentric(MARS, t)
    t0 = time.time()
    chain = naff_vector(X, t, 24)
    assert time.time() - t0 < 20
    ref = [c for c in geocentric_chain(MARS, EARTH, 40, 40, t_ref=t[0]) if c.omega != 0]
    ref.sort(key=lambda c: -c.radius)
    fe = re = 0.0
    for e in ref[:10]:
        Om = _canon(e)
        c = min(chain, key=lambda q: np.linalg.norm(_canon(q) - Om))
        fe = max(fe, np.linalg.norm(_canon(c) - Om) / np.linalg.norm(Om))
        re = max(re, abs(c.radius / e.radius - 1))
    # measured at 10-day sampling, K=24: fe = 2.0e-10, re = 8.9e-8, reconstruction 1.9e-8 AU
    assert fe < 1e-8, fe
    assert re < 1e-6, re
    # fixed arm: sum of the analytic (omega = 0) arms
    arms = sum(c.radius * c.u * np.cos(c.phase) + c.radius * c.v * np.sin(c.phase)
               for c in geocentric_chain(MARS, EARTH, 40, 40, t_ref=t[0]) if c.omega == 0)
    fixed = [c for c in chain if c.omega == 0]
    assert len(fixed) == 1 and fixed[0].radius == pytest.approx(np.linalg.norm(arms), rel=1e-6)
    # reconstruction by evaluate_chain of the NAFF chain
    err = np.max(np.linalg.norm(evaluate_chain(chain, t) - X, axis=1))
    assert err < 1e-6, err        # AU, max over the record (largest at the ends where w -> 0)
