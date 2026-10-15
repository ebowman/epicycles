"""Export golden test data for the TypeScript port (web/test/golden/*.json).

Run: uv run python experiments/export_golden.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv, jvp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from epicycles.bodies import EARTH, MARS, PHAETHON  # noqa: E402
from epicycles.expansion import (chain_arms, evaluate_chain, geocentric_chain,  # noqa: E402
                                 kepler_fourier)
from epicycles.export import chain_to_json  # noqa: E402
from epicycles.kepler import position, solve_kepler  # noqa: E402

OUT = ROOT / "web" / "test" / "golden"
OUT.mkdir(parents=True, exist_ok=True)


def dump(name, obj):
    (OUT / name).write_text(json.dumps(obj, separators=(",", ":")))


# --- solve_kepler
kep = []
for e in (0.0, 0.0167, 0.0934, 0.5, 0.8, 0.89, 0.99):
    for M in (-20.0, -3.1, -1e-3, 0.0, 1e-6, 0.7, 2.0, 3.14159, 6.0, 100.5):
        kep.append({"M": M, "e": e, "E": solve_kepler(M, e)})
dump("kepler.json", kep)

# --- bessel
bes = []
for n in (0, 1, 2, 5, 10, 30, 60, 100, 200, 300, 450, 600, 601):
    xs = {0.0, 0.1, 1.0, n * 0.5, n * 0.9, float(n), n * 1.1 if n else 2.0, n + 3.0, 7.5}
    if n:
        xs |= {n * 0.89, n * 0.99}
    for x in sorted(xs):
        if x > 650:
            continue
        bes.append({"n": n, "x": x, "J": float(jv(n, x)), "Jp": float(jvp(n, x))})
dump("bessel.json", bes)

# --- kepler_fourier
kf = []
for e, Ns in ((0.0, (0, 1, 5, 60)), (1e-10, (1, 5, 60)), (0.0934, (5, 20, 60)),
              (0.5, (5, 60)), (0.89, (20, 60, 300))):
    for N in Ns:
        kf.append({"e": e, "N": N, "c": kepler_fourier(e, N).real.tolist()})
dump("kepler_fourier.json", kf)

# --- positions
pos = []
for el in (MARS, EARTH, PHAETHON):
    for dt in (-5000.0, -123.4, 0.0, 17.25, 400.0, 3000.0, 9000.5):
        t = el.epoch + dt
        pos.append({"body": el.name.upper(), "t": t, "r": position(el, t).tolist()})
dump("positions.json", pos)

# --- geocentric chains
gc = []
for body, nb, ne in ((MARS, 5, 5), (MARS, 20, 20), (PHAETHON, 60, 60)):
    chain = geocentric_chain(body, EARTH, nb, ne)
    for dt in (-4000.0, -10.5, 0.0, 33.3, 700.0, 5000.0):
        t = body.epoch + dt
        gc.append({"body": body.name.upper(), "Nb": nb, "Ne": ne, "t": t,
                   "r": evaluate_chain(chain, t).tolist(),
                   "truth": (position(body, t) - position(EARTH, t)).tolist(),
                   "arms_last": chain_arms(chain, t)[-1].tolist(),
                   "arms_mid": chain_arms(chain, t)[len(chain) // 2].tolist()})
dump("geocentric.json", gc)

# --- full chain (field by field)
chain = geocentric_chain(MARS, EARTH, 5, 5)
dump("mars_N5_chain.json", chain_to_json(chain, {"body": "Mars", "earth": "Earth", "N_body": 5,
                                                  "N_earth": 5, "t_ref": MARS.epoch}))
print("wrote golden files to", OUT)
