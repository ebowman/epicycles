"""Generating script for Table I (tab:nneeded) and Table IV (tab:ptolemy) of the paper.

Run: uv run python experiments/tables.py     (fast, a few seconds)
Writes experiments/out/nneeded.json and experiments/out/ptolemy_truncations.json.

Table I.  For each e: smallest N with sampled max_M |z - z_N|/a < tol (3001-point grid on
[0, pi]; the error is symmetric under M -> -M), tol in {1e-6, 1e-12}, and the smallest N with
epicycles.bounds.tail_bound(e, N) < tol.  z is the perifocal position x + iy, focus at origin.
Table IV.  Planar (coplanar) geocentric Mars, J2000 mean elements, 60 yr at 0.5 d steps, exact
two-body reference; longitude error of harmonic truncations A-G and the retrograde-interval count.
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from epicycles.bodies import EARTH, MARS  # noqa: E402
from epicycles.bounds import rho, tail_bound  # noqa: E402
from epicycles.expansion import evaluate_chain, geocentric_chain, kepler_fourier  # noqa: E402
from epicycles.frames import ecliptic_lonlat, geocentric  # noqa: E402
from epicycles.kepler import solve_kepler  # noqa: E402

OUT = Path(__file__).parent / "out"
ES = [(0.0167, "Earth"), (0.0934, "Mars"), (0.25, ""), (0.6, ""), (0.89, "")]
TOLS = [1e-6, 1e-12]
NGRID = 3001


def _first_below(errs, tol):
    idx = np.nonzero(errs < tol)[0]
    return int(idx[0]) if len(idx) else None


def n_needed(e, nmax, ngrid=NGRID):
    """Return (observed sup-error array indexed by N=0..nmax)."""
    M = np.linspace(0.0, np.pi, ngrid)
    E = solve_kepler(M, e)
    z = (np.cos(E) - e) + 1j * np.sqrt(1 - e * e) * np.sin(E)
    c = kepler_fourier(e, nmax)                       # index k + nmax
    zN = c[nmax] * np.ones_like(M, dtype=complex)
    errs = [float(np.max(np.abs(z - zN)))]
    for k in range(1, nmax + 1):
        zN = zN + c[nmax + k] * np.exp(1j * k * M) + c[nmax - k] * np.exp(-1j * k * M)
        errs.append(float(np.max(np.abs(z - zN))))
    return np.array(errs)


def table_nneeded():
    rows = []
    for e, label in ES:
        errs = n_needed(e, 1500 if e > 0.8 else 400)
        row = {"e": e, "label": label, "rho": round(rho(e), 4), "observed": {}, "bound": {}}
        for tol in TOLS:
            row["observed"][f"{tol:g}"] = _first_below(errs, tol)
            bn = 0
            while tail_bound(e, bn) >= tol:
                bn += 1
            row["bound"][f"{tol:g}"] = bn
        rows.append(row)
    return {"grid_points": NGRID, "tols": TOLS, "rows": rows}


def _planar(el):
    return replace(el, i=0.0, Omega=0.0, omega=el.omega + el.Omega)


def _angdiff_deg(a, b):
    return np.degrees(np.abs((a - b + np.pi) % (2 * np.pi) - np.pi))


MODELS = [  # name, description, body k-set, earth k-set, Nmax
    ("A", "k=1, both bodies", {1}, {1}),
    ("B", "k in {0,1}, both", {0, 1}, {0, 1}),
    ("C", "Mars k in {0,1,2}, Earth k=1", {0, 1, 2}, {1}),
    ("D", "k in {0,1,2}, both", {0, 1, 2}, {0, 1, 2}),
    ("E", "|k|<=2, both", set(range(-2, 3)), set(range(-2, 3))),
    ("F", "|k|<=5, both", set(range(-5, 6)), set(range(-5, 6))),
    ("G", "|k|<=10, both", set(range(-10, 11)), set(range(-10, 11))),
]


def table_ptolemy():
    mp, ep = _planar(MARS), _planar(EARTH)
    t = mp.epoch + np.arange(0, 60 * 365.25, 0.5)
    lon_t, _ = ecliptic_lonlat(geocentric(mp, t, ep))

    def n_retro(lon):
        """Number of maximal runs of decreasing unwrapped longitude."""
        r = np.diff(np.unwrap(lon)) < 0
        return int(np.sum(np.diff(r.astype(int)) == 1) + int(r[0]))

    rows = []
    for name, desc, kb, ke in MODELS:
        nmax = max(max(abs(k) for k in kb), max(abs(k) for k in ke))
        chain = geocentric_chain(mp, ep, nmax, nmax, t_ref=mp.epoch)
        sel = [c for c in chain if (c.source == "body" and c.k in kb) or (c.source == "earth" and c.k in ke)]
        lon, _ = ecliptic_lonlat(evaluate_chain(sel, t))
        d = _angdiff_deg(lon, lon_t)
        rows.append({"model": name, "terms": desc, "max_deg": float(d.max()),
                     "rms_deg": float(np.sqrt((d ** 2).mean())), "n_retro_intervals": n_retro(lon)})
    return {"span_years": 60, "dt_days": 0.5, "n_retro_true": n_retro(lon_t), "rows": rows}


def main():
    OUT.mkdir(exist_ok=True)
    t1, t4 = table_nneeded(), table_ptolemy()
    (OUT / "nneeded.json").write_text(json.dumps(t1, indent=1) + "\n")
    (OUT / "ptolemy_truncations.json").write_text(json.dumps(t4, indent=1) + "\n")
    for r in t1["rows"]:
        print(r["e"], r["rho"], *r["observed"].values(), *r["bound"].values())
    print("true retro intervals:", t4["n_retro_true"])
    for r in t4["rows"]:
        print(r["model"], f"{r['max_deg']:.4g} {r['rms_deg']:.4g}", r["n_retro_intervals"])


if __name__ == "__main__":
    main()
