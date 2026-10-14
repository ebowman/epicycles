"""JSON export of epicycle chains (schema: docs/epicycle-schema.md)."""
from __future__ import annotations

import json
from pathlib import Path

from .expansion import Epicycle

SCHEMA_VERSION = 1


def chain_to_json(chain: list[Epicycle], meta: dict) -> dict:
    """Build the schema-v1 dict.  meta needs body, earth, N_body, N_earth and
    optionally t_ref (JD; default: chain[0].t_ref).  All circles are
    re-referenced to meta t_ref."""
    if not chain and "t_ref" not in meta:
        raise ValueError("t_ref required for an empty chain")
    t_ref = float(meta.get("t_ref", chain[0].t_ref if chain else 0.0))
    m = {"body": meta.get("body"), "earth": meta.get("earth"),
         "N_body": meta.get("N_body"), "N_earth": meta.get("N_earth"),
         "units": {"length": "AU", "time": "day", "angle": "rad"}, "t_ref": t_ref}
    eps = []
    for c in chain:
        c = c.rereferenced(t_ref)
        eps.append({"radius": float(c.radius), "omega": float(c.omega), "phase": float(c.phase),
                    "normal": [float(x) for x in c.normal],
                    "u": [float(x) for x in c.u], "v": [float(x) for x in c.v],
                    "source": c.source, "k": c.k})
    return {"schema_version": SCHEMA_VERSION, "meta": m, "epicycles": eps}


def write_chain(path, chain: list[Epicycle], meta: dict) -> None:
    """Write chain_to_json(chain, meta) to path (parents created)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(chain_to_json(chain, meta), indent=1) + "\n")


def chain_from_json(d: dict) -> list[Epicycle]:
    """Inverse of chain_to_json (for round-trip tests / Python consumers)."""
    import numpy as np
    t_ref = d["meta"]["t_ref"]
    return [Epicycle(e["radius"], e["omega"], e["phase"], np.array(e["normal"]),
                     np.array(e["u"]), np.array(e["v"]), t_ref, e["source"], e["k"])
            for e in d["epicycles"]]
