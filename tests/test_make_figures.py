import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))
import make_figures as mf  # noqa: E402

from epicycles.bodies import EARTH, MARS  # noqa: E402
from epicycles.frames import geocentric  # noqa: E402


def eval_json(d, t):
    t = np.atleast_1d(np.asarray(t, float))
    r = np.zeros((len(t), 3))
    for e in d["epicycles"]:
        th = e["phase"] + e["omega"] * (t - d["meta"]["t_ref"])
        r += e["radius"] * (np.cos(th)[:, None] * np.array(e["u"]) + np.sin(th)[:, None] * np.array(e["v"]))
    return r


def test_export_data(tmp_path):
    index = mf.export_data(tmp_path)
    idx = json.loads((tmp_path / "index.json").read_text())["files"]
    assert idx == index
    names = {e["file"] for e in idx}
    for b in ("mars", "ceres", "phaethon"):
        for N in (1, 5, 20, 60):
            assert f"{b}_N{N}.json" in names
    assert "perturbed_ceres.json" in names
    assert {e["kind"] for e in idx} == {"analytic", "naff"}
    for e in idx:
        d = json.loads((tmp_path / e["file"]).read_text())
        assert d["schema_version"] == 1
        assert d["meta"]["units"] == {"length": "AU", "time": "day", "angle": "rad"}
        assert "t_ref" in d["meta"] and len(d["epicycles"]) == e["n_circles"]
        for c in d["epicycles"][:5]:
            assert set(c) >= {"radius", "omega", "phase", "normal", "u", "v", "source", "k"}
            assert c["source"] in ("body", "earth", "naff")
        radii = [c["radius"] for c in d["epicycles"]]
        assert radii == sorted(radii, reverse=True)
    pc = json.loads((tmp_path / "perturbed_ceres.json").read_text())
    assert "study" in pc and {c["source"] for c in pc["epicycles"]} == {"naff"}

    mars = json.loads((tmp_path / "mars_N60.json").read_text())
    t = MARS.epoch + np.array([0.0, 123.4, 2000.0, 9000.5])
    assert np.max(np.abs(eval_json(mars, t) - geocentric(MARS, t, EARTH))) < 1e-9


def test_smoke_fig(tmp_path, monkeypatch):
    monkeypatch.setattr(mf, "FIG_DIR", tmp_path)
    mf.style()
    rho = mf.fig_coeff_decay()
    assert (tmp_path / "fig_coeff_decay.pdf").stat().st_size > 1000
    assert abs(rho[0.0934] - 0.1267) < 1e-3
