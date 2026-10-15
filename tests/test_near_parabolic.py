"""Tests for experiments/near_parabolic.py (review item 1): limiting coefficients, N_actual vs e."""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv, jvp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))
import near_parabolic as npb  # noqa: E402

from epicycles.expansion import kepler_fourier  # noqa: E402


def test_limit_coefficient_is_jprime_over_k():
    k = np.array([1.0, 2.0, 5.0, 50.0, 1000.0, 1e5])
    cp, cm = npb.coeff_pair(1.0, k)
    assert np.allclose(cp, jvp(k, k) / k, rtol=1e-9, atol=0)
    assert np.allclose(cp, cm)
    assert np.allclose(cp, (jv(k, k) - jv(k + 1, k)) / k, rtol=1e-9, atol=0)


def test_limit_coefficient_asymptotics():
    k = np.array([1e4, 1e6])
    cp, _ = npb.coeff_pair(1.0, k)
    assert np.allclose(cp * k ** (5 / 3), npb.C_LIM, rtol=2e-3)
    assert abs(npb.C_LIM - 0.41085) < 1e-5


def test_coefficients_approach_limit_as_e_to_1():
    k = np.arange(1, 6)
    c = kepler_fourier(1 - 1e-14, 5)
    lim, _ = npb.coeff_pair(1.0, k)
    assert np.allclose(c[6:], lim, atol=1e-6) and np.allclose(c[:5][::-1], lim, atol=1e-6)


def test_tail_at_M0_equals_sampled_sup_and_signs():
    e, N = 0.9, 40
    cross, minc, firstneg, *_ = npb.stream_tail(e, [1e-2], kmax=200000, grid=False)
    assert minc >= 0 and firstneg is None            # all coefficients nonnegative -> sup at M=0
    mx, Mmx, T0 = npb.sampled_tail_max(e, N, nM=800)
    assert mx == np.float64(mx) and abs(mx - T0) < 1e-12 and Mmx < 1e-6


def test_n_actual_monotone_in_e_and_below_cert():
    es = [0.5, 0.89, 0.99, 0.999, 0.9999]
    na, nc = [], []
    for e in es:
        cross, *_ = npb.stream_tail(e, [1e-2], kmax=3_000_000, grid=False)
        na.append(cross[1e-2]); nc.append(npb.n_cert(e, 1e-2))
    assert all(a is not None for a in na)
    assert all(x < y for x, y in zip(na, na[1:]))
    assert all(a <= c for a, c in zip(na, nc))
    cross1, *_ = npb.stream_tail(1.0, [1e-2], kmax=3_000_000, grid=False)
    assert na[-1] <= cross1[1e-2] <= 1368 + 1       # finite limit near 1367
