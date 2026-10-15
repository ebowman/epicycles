import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))
import tables as tb  # noqa: E402


def test_table_nneeded_matches_paper():
    rows = tb.table_nneeded()["rows"]
    obs = [(r["observed"]["1e-06"], r["observed"]["1e-12"]) for r in rows]
    bnd = [(r["bound"]["1e-06"], r["bound"]["1e-12"]) for r in rows]
    assert obs == [(4, 7), (6, 12), (10, 22), (33, 75), (230, 575)]
    assert bnd == [(4, 8), (7, 13), (12, 24), (43, 87), (350, 712)]


def test_table_ptolemy_matches_paper():
    t = tb.table_ptolemy()
    assert t["n_retro_true"] == 28
    r = {x["model"]: x for x in t["rows"]}
    assert all(x["n_retro_intervals"] == 28 for x in r.values())
    assert abs(r["A"]["max_deg"] - 33.2) < 0.05 and abs(r["B"]["max_deg"] - 9.40) < 0.005
    assert abs(r["C"]["max_deg"] - 5.05) < 0.005 and abs(r["D"]["max_deg"] - 0.557) < 0.0005
    assert abs(r["E"]["max_deg"] - 0.728) < 0.0005 and r["F"]["max_deg"] < 7e-4
    assert r["G"]["max_deg"] < 1e-4
