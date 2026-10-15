"""Convergence of NAFF epicycle extraction on an N-body-perturbed geocentric Ceres.

Run: uv run python experiments/perturbed_ceres.py
Writes experiments/out/perturbed_ceres.json and perturbed_ceres_errors.csv.
"""
import csv
import json
from pathlib import Path

import numpy as np
from dataclasses import replace
from scipy.optimize import least_squares

from epicycles.bodies import CERES, EARTH, JUPITER, SATURN
from epicycles.expansion import evaluate_chain, geocentric_chain
from epicycles.export import chain_to_json
from epicycles.naff import hann_window, naff_vector
from epicycles.kepler import MU_SUN, position
from epicycles.nbody import MASS_RATIOS, geocentric, integrate

OUT = Path(__file__).parent / "out"
YEARS, DT, WINDOW_P = 200, 5.0, 4
KS = [5, 10, 20, 50, 100, 200]
KEPLER_N = [10, 20, 40]


def errs(truth, rec):
    d = np.linalg.norm(truth - rec, axis=1)
    n = len(d)
    i = d[int(0.1 * n):int(0.9 * n)]
    return {"full_max": float(d.max()), "full_rms": float(np.sqrt(np.mean(d**2))),
            "interior_max": float(i.max()), "interior_rms": float(np.sqrt(np.mean(i**2)))}


def fit_kepler(el0, t, g):
    """Least-squares 6-element Kepler orbit for Ceres: min |Ceres_kep - EARTH_kep - g|
    over the window (every 5th sample), starting from el0."""
    idx = slice(None, None, 5)
    ref = position(EARTH, t[idx])

    def mk(x):
        return replace(el0, a=x[0], e=x[1], i=x[2], Omega=x[3], omega=x[4], M0=x[5])

    def res(x):
        return ((position(mk(x), t[idx]) - ref) - g[idx]).ravel()

    x0 = [el0.a, el0.e, el0.i, el0.Omega, el0.omega, el0.M0]
    sc = [1.0, 0.1, 0.1, 0.1, 0.1, 0.1]
    r = least_squares(res, x0, x_scale=sc, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return mk(r.x)


def fit_kepler_helio(el0, t, h, weight=None):
    """6-element Kepler fit to HELIOCENTRIC Ceres h (every 5th sample).  ``weight``
    (per-sample, same length as t) multiplies the residuals by sqrt(weight) -- pass the
    NAFF Hann^p window to match NAFF's weighting; None = unweighted."""
    idx = slice(None, None, 5)
    sw = 1.0 if weight is None else np.sqrt(weight[idx] / weight.max())[:, None]

    def mk(x):
        return replace(el0, a=x[0], e=x[1], i=x[2], Omega=x[3], omega=x[4], M0=x[5])

    def res(x):
        return (sw * (position(mk(x), t[idx]) - h[idx])).ravel()

    x0 = [el0.a, el0.e, el0.i, el0.Omega, el0.omega, el0.M0]
    sc = [1.0, 0.1, 0.1, 0.1, 0.1, 0.1]
    r = least_squares(res, x0, x_scale=sc, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return mk(r.x)


def main():
    bodies = {"Earth": EARTH, "Jupiter": JUPITER, "Saturn": SATURN, "Ceres": CERES}
    t0 = CERES.epoch
    tr = integrate(bodies, t0, YEARS * 365.25, DT)
    t, g = tr["t"], geocentric(tr, "Ceres", "Earth")
    print(f"N={len(t)} samples, {YEARS} yr, dt={DT} d, epoch JD {t0}")
    res, best = {}, None
    rows = []
    for K in KS:
        ch = naff_vector(g, t, K, window_p=WINDOW_P)
        e = errs(g, evaluate_chain(ch, t))
        res[K] = {"n_circles": len(ch), **e}
        rows.append(("naff", K, len(ch), e))
        best = ch
    kep = {}
    kep_other = {}
    for N in KEPLER_N:
        e = errs(g, evaluate_chain(geocentric_chain(CERES, EARTH, N, N, t_ref=t0), t))
        kep[N] = e
        rows.append(("kepler_chain", N, None, e))
    # (b) mean-element Kepler: a from the dominant NAFF line of heliocentric Ceres
    hel = naff_vector(tr["Ceres"], t, 3, window_p=WINDOW_P)
    dom = max((c for c in hel if c.omega != 0), key=lambda c: c.radius)
    n_mean = abs(dom.omega)
    mean_el = replace(CERES, a=(MU_SUN / n_mean**2) ** (1 / 3))
    # (c) best-fit 6-element Kepler orbit
    fit_el = fit_kepler(mean_el, t, g)
    for name, el in (("kepler_mean_a", mean_el), ("kepler_fit6", fit_el)):
        # exact Kepler (N -> infinity) geocentric position
        e = errs(g, position(el, t) - position(EARTH, t))
        kep_other[name] = {**e, "a": el.a, "e": el.e}
        rows.append((name, 0, None, e))
    # ---- fair baseline: heliocentric Ceres h = r_Ceres - r_Sun (integrated), integrated
    # Earth E (heliocentric); every reconstruction is differenced with the SAME integrated E,
    # so Earth-model error is excluded: err = |(recon_h - E) - g|.
    h, E = tr["Ceres"], tr["Earth"]
    hel_k, hel_series = {}, {}
    for K in KS:
        rec = evaluate_chain(naff_vector(h, t, K, window_p=WINDOW_P), t)
        e = errs(g, rec - E)
        hel_k[K] = e
        rows.append(("naff_helio", K, 2 * K + 1, e))
        hel_series[f"naff_helio_K{K}"] = np.linalg.norm(rec - E - g, axis=1)
    w = hann_window(len(t), WINDOW_P)
    hel_fit = {}
    for name, wt in (("kepler_helio_unweighted", None), ("kepler_helio_hann", w)):
        el = fit_kepler_helio(mean_el, t, h, wt)
        rec = position(el, t) - E
        e = errs(g, rec)
        hel_fit[name] = {**e, "a": el.a, "e": el.e,
                         "helio_full_max": float(np.linalg.norm(position(el, t) - h, axis=1).max())}
        rows.append((name, 0, None, e))
        hel_series[name] = np.linalg.norm(rec - g, axis=1)
    OUT.mkdir(exist_ok=True)
    np.savez_compressed(OUT / "perturbed_ceres_residuals.npz", t=t, **hel_series)
    # ---- tolerance check of the actual 200-yr integration (rtol 1e-12/atol 1e-14 vs tighter)
    tr2 = integrate(bodies, t0, YEARS * 365.25, DT, rtol=1e-13, atol=1e-16)
    g2 = geocentric(tr2, "Ceres", "Earth")
    tol = {"baseline": {"method": "DOP853", "rtol": 1e-12, "atol": 1e-14},
           "tight": {"method": "DOP853", "rtol": 1e-13, "atol": 1e-16},
           "max_dg_AU": float(np.linalg.norm(g2 - g, axis=1).max()),
           "max_dh_ceres_AU": float(np.linalg.norm(tr2["Ceres"] - h, axis=1).max()),
           "max_dEarth_AU": float(np.linalg.norm(tr2["Earth"] - E, axis=1).max()),
           "note": "max over the 200-yr sampled grid of |g_tight - g_baseline|"}
    (OUT / "perturbed_ceres_tolerance.json").write_text(json.dumps(tol, indent=1) + "\n")
    print("tolerance check:", tol)
    print(f"{'method':14s}{'K/N':>5s}{'circ':>6s}{'full_max':>11s}{'full_rms':>11s}{'int_max':>11s}{'int_rms':>11s}")
    for m, K, nc, e in rows:
        print(f"{m:14s}{K:5d}{(nc or 0):6d}{e['full_max']:11.2e}{e['full_rms']:11.2e}"
              f"{e['interior_max']:11.2e}{e['interior_rms']:11.2e}")
    meta = {"body": "Ceres (N-body perturbed)", "earth": "Earth-Moon barycentre (N-body)",
            "N_body": None, "N_earth": None, "t_ref": float(t[0])}
    doc = chain_to_json(best, meta)
    doc["study"] = {
        "masses_Msun": {"Sun": 1.0, **{k: v for k, v in MASS_RATIOS.items()}, "Ceres": 0.0},
        "epoch_jd": t0, "span_years": YEARS, "dt_days": DT, "n_samples": len(t),
        "window": f"Hann^p, p={WINDOW_P} (sin^{2 * WINDOW_P})", "K_max": KS[-1], "chain_K": KS[-1],
        "interior": "middle 80% of samples", "units": "AU",
        "naff_errors": {str(k): v for k, v in res.items()},
        "kepler_osculating_chain_errors": {str(k): v for k, v in kep.items()},
        "kepler_other": kep_other,
        "heliocentric_baseline": {
            "naff_helio_errors": {str(k): v for k, v in hel_k.items()},
            "kepler_helio_fits": hel_fit,
            "note": "NAFF (Hann^p) on heliocentric Ceres and 6-element Kepler fits to heliocentric "
                    "Ceres (unweighted; weighted by the NAFF Hann^p window, residual weights "
                    "sqrt(w)); geocentric error = |(recon_h - E_int) - g| with the SAME integrated "
                    "Earth E_int subtracted, so Earth-model error is excluded."},
        "kepler_note": "osculating-element Kepler error is dominated by the osculating-vs-mean "
                       "element offset and is NOT a fair Kepler-vs-epicycle comparison; "
                       "kepler_mean_a uses a from the dominant NAFF line of heliocentric Ceres; "
                       "kepler_fit6 is a least-squares 6-element ellipse (exact Kepler, Earth Kepler).",
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "perturbed_ceres.json").write_text(json.dumps(doc, indent=1) + "\n")
    with open(OUT / "perturbed_ceres_errors.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "K_or_N", "n_circles", "full_max", "full_rms", "interior_max", "interior_rms"])
        for m, K, nc, e in rows:
            w.writerow([m, K, nc, e["full_max"], e["full_rms"], e["interior_max"], e["interior_rms"]])


if __name__ == "__main__":
    main()
