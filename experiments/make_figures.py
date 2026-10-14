"""Single deterministic entry point: publication figures + SPA data exports.

Run:  uv run python experiments/make_figures.py
Writes figures/*.pdf, figures/captions.md and web/public/data/*.json.
Reads (does not re-run) experiments/out/perturbed_ceres.json -- regenerate that
first with `uv run python experiments/perturbed_ceres.py` (~3 min) if missing.
"""
from __future__ import annotations

import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from epicycles.bodies import CERES, EARTH, EROS, MARS, PHAETHON  # noqa: E402
from epicycles.bounds import geocentric_tail_bound, rho  # noqa: E402
from epicycles.expansion import (chain_arms, evaluate_chain, geocentric_chain,  # noqa: E402
                                 kepler_fourier)
from epicycles.export import chain_to_json  # noqa: E402
from epicycles.frames import ecliptic_lonlat, geocentric  # noqa: E402
from epicycles.kepler import mean_motion  # noqa: E402

FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "web" / "public" / "data"
PERTURBED_JSON = ROOT / "experiments" / "out" / "perturbed_ceres.json"

# Okabe-Ito colorblind-safe palette
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "yellow": "#F0E442", "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7"}
BODY_COLOR = {"MARS": OI["verm"], "CERES": OI["blue"], "PHAETHON": OI["green"]}
W1 = 3.4  # single-column width [in]


def style():
    plt.rcParams.update({
        "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.5,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6,
        "lines.linewidth": 1.1, "legend.frameon": False, "figure.dpi": 150,
        "axes.prop_cycle": plt.cycler(color=[OI[k] for k in
                                             ("blue", "verm", "green", "orange", "purple", "sky")]),
        "pdf.fonttype": 42,
    })


def save(fig, name):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / name, bbox_inches="tight", metadata={"CreationDate": None, "Creator": None})
    plt.close(fig)


# ---------------------------------------------------------------- (a)
def fig_coeff_decay() -> dict:
    N = 150
    es = [0.0934, 0.25, 0.6, 0.89]
    cols = [OI["blue"], OI["green"], OI["orange"], OI["verm"]]
    fig, ax = plt.subplots(figsize=(W1, 3.4))
    k = np.arange(1, N + 1)
    info = {}
    mk = ["o", "s", "^", "D"]
    for e, col, m in zip(es, cols, mk):
        c = kepler_fourier(e, N)
        pos, neg = np.abs(c[N + 1:]), np.abs(c[:N][::-1])
        ax.semilogy(k, pos, m, ms=2.2, mfc="none", mew=0.6, color=col, label=f"$e={e}$, $|c_k|$")
        ax.semilogy(k, neg, m, ms=2.2, mew=0.6, color=col, label=f"$e={e}$, $|c_{{-k}}|$")
        ax.semilogy(k, rho(e) ** k, ls=(0, (4, 1.5 + 1.5 * len(info))), lw=0.8, color=col)
        info[e] = rho(e)
    ax.set_ylim(1e-16, 3)
    ax.set_xlim(0, N)
    ax.set_xlabel("harmonic index $k$")
    ax.set_ylabel(r"coefficient magnitude $|c_{\pm k}|$ (units of $a$)")
    ax.plot([], [], "--", color="gray", lw=0.8, label=r"$\rho(e)^k$ (Kapteyn)")
    ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17), handlelength=1.2, columnspacing=0.8)
    save(fig, "fig_coeff_decay.pdf")
    return info


# ---------------------------------------------------------------- (b)
def chain_error(body, N, t, truth):
    ch = geocentric_chain(body, EARTH, N, N, t_ref=float(t[0]))
    return float(np.max(np.linalg.norm(evaluate_chain(ch, t) - truth, axis=1)))


def fig_geo_convergence() -> dict:
    t = 2461200.5 + np.arange(0.0, 50 * 365.25, 1.0)
    Ns = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200,
          210, 220, 230, 240, 250, 260, 280, 300, 325, 350, 375, 400, 450, 500, 600]
    fig, ax = plt.subplots(figsize=(W1, 2.7))
    res = {}
    for body in (MARS, CERES, PHAETHON):
        truth = geocentric(body, t)
        err = [chain_error(body, N, t, truth) for N in Ns]
        bnd = [geocentric_tail_bound(body, EARTH, N, N) for N in Ns]
        col = BODY_COLOR[body.name.upper()]
        lab = f"{body.name} ($e={body.e:.3f}$)"
        ax.loglog(Ns, np.maximum(err, 1e-17), "o-", ms=2.2, color=col, label=lab)
        ax.loglog(Ns, np.maximum(bnd, 1e-17), "--", lw=0.8, color=col)
        res[body.name.upper()] = dict(
            Ns=Ns, err=err, bound=bnd,
            N_1e6=next((n for n, x in zip(Ns, err) if x < 1e-6), None),
            N_1e9=next((n for n, x in zip(Ns, err) if x < 1e-9), None),
            floor=min(err))
    floor = 2e-13
    ax.axhline(floor, color="0.4", lw=0.7, ls=":")
    ax.text(25, floor * 1.5, "roundoff floor (~2e-13 AU)", fontsize=6, color="0.3")
    ax.plot([], [], "--", color="gray", lw=0.8, label="tail bound (Part 4)")
    ax.set_ylim(3e-14, 5)
    ax.set_xlabel(r"truncation order $N$ ($N_{\rm body}=N_\oplus=N$)")
    ax.set_ylabel("max geocentric error [AU] (50 yr)")
    ax.legend(loc="lower left", bbox_to_anchor=(0.42, 0.12), handlelength=1.4)
    save(fig, "fig_geo_convergence.pdf")
    return res


# ---------------------------------------------------------------- (c)
def _planar(el):
    return replace(el, i=0.0, Omega=0.0, omega=el.omega + el.Omega)


def _angdiff_deg(a, b):
    return np.degrees(np.abs((a - b + np.pi) % (2 * np.pi) - np.pi))


def fig_mars_retrograde() -> dict:
    # Ptolemaic order = derivation Part 3 / 4.3(d) model A: only k=+1 of Mars (deferent)
    # and k=+1 of Earth (epicycle), 3D here; the derivation's planar 60-yr statistic
    # is recomputed below so the caption numbers agree with the theory note.
    mp, ep_ = _planar(MARS), _planar(EARTH)
    tp = mp.epoch + np.arange(0, 60 * 365.25, 0.5)
    lon_t, _ = ecliptic_lonlat(geocentric(mp, tp, ep_))
    pt_p = [c for c in geocentric_chain(mp, ep_, 1, 1, t_ref=mp.epoch) if c.k == 1]
    lon_p, _ = ecliptic_lonlat(evaluate_chain(pt_p, tp))
    d = _angdiff_deg(lon_p, lon_t)
    planar_stats = {"max_deg": float(d.max()), "rms_deg": float(np.sqrt((d**2).mean()))}

    t0 = 2459000.5  # 2020-05-31; spans the 2020-10 and 2022-12 oppositions
    t = t0 + np.arange(0, 2.75 * 365.25, 0.5)
    truth = geocentric(MARS, t)
    pt = [c for c in geocentric_chain(MARS, EARTH, 1, 1, t_ref=t0) if c.k == 1]
    n5 = geocentric_chain(MARS, EARTH, 5, 5, t_ref=t0)
    curves = {"true (Kepler)": truth, "Ptolemaic ($k=1$ deferent + epicycle)": evaluate_chain(pt, t),
              "$N=5$": evaluate_chain(n5, t)}
    ll = {k: ecliptic_lonlat(v) for k, v in curves.items()}
    lon_u = {k: np.degrees(np.unwrap(v[0])) for k, v in ll.items()}
    styles = {"true (Kepler)": dict(color="k", lw=1.6),
              "Ptolemaic ($k=1$ deferent + epicycle)": dict(color=OI["verm"], lw=1.0, ls="--"),
              "$N=5$": dict(color=OI["sky"], lw=1.0)}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(W1, 4.6), height_ratios=[1.2, 1])
    for k, (lon, lat) in ll.items():
        x = np.degrees(lon)
        brk = np.where(np.abs(np.diff(x)) > 180)[0] + 1
        xs, ys = np.insert(x, brk, np.nan), np.insert(np.degrees(lat), brk, np.nan)
        a1.plot(xs, ys, label=k, **styles[k])
    a1.set_xlabel("ecliptic longitude [deg]")
    a1.set_ylabel("ecliptic latitude [deg]")
    a1.legend(loc="best", handlelength=1.6)
    yrs = (t - t0) / 365.25
    dlon = np.diff(lon_u["true (Kepler)"])
    retro = dlon < 0
    for k in curves:
        a2.plot(yrs, lon_u[k] - 0.5240 * (t - t0), label=k, **styles[k])
    ymin, ymax = a2.get_ylim()
    a2.fill_between(yrs[1:], ymin, ymax, where=retro, color="0.85", step="mid", lw=0, zorder=0)
    a2.set_ylim(ymin, ymax)
    a2.set_xlim(yrs[0], yrs[-1])
    a2.set_xlabel("time since 2020-05-31 [yr]")
    a2.set_ylabel("longitude minus mean motion [deg]")
    fig.tight_layout()
    save(fig, "fig_mars_retrograde.pdf")
    e_pt = _angdiff_deg(ll["Ptolemaic ($k=1$ deferent + epicycle)"][0], ll["true (Kepler)"][0])
    e_5 = _angdiff_deg(ll["$N=5$"][0], ll["true (Kepler)"][0])
    return {"planar_60yr": planar_stats, "window_ptol_max_deg": float(e_pt.max()),
            "window_N5_max_deg": float(e_5.max()),
            "n_retro_intervals_window": int(np.sum(np.diff(retro.astype(int)) == 1) + int(retro[0]))}


# ---------------------------------------------------------------- (d)
def fig_chain_3d() -> dict:
    N = 3
    t_ref = EROS.epoch
    chain = geocentric_chain(EROS, EARTH, N, N, t_ref=t_ref)
    syn = 2 * np.pi / abs(mean_motion(EROS) - mean_motion(EARTH))
    ts = t_ref + np.linspace(0, syn, 500)
    trail = evaluate_chain(chain, ts)
    truth = geocentric(EROS, ts)
    t_snap = t_ref + 0.62 * syn
    arms = chain_arms(chain, t_snap)
    fig = plt.figure(figsize=(W1, 3.5))
    ax = fig.add_subplot(projection="3d")
    ax.plot(*truth.T, color="0.7", lw=2.4, label="true path", zorder=1)
    ax.plot(*trail.T, color=OI["blue"], lw=0.9, label=f"chain tip, $N={N}$", zorder=2)
    th = np.linspace(0, 2 * np.pi, 120)
    for j, c in enumerate(chain):
        col = OI["verm"] if c.source == "body" else OI["green"]
        ctr = arms[j]
        ring = ctr + c.radius * (np.cos(th)[:, None] * c.u + np.sin(th)[:, None] * c.v)
        ax.plot(*ring.T, color=col, lw=0.5, alpha=0.8)
        ax.plot(*arms[j:j + 2].T, color=col, lw=1.0)
    ax.plot([], [], color=OI["verm"], label="body circles")
    ax.plot([], [], color=OI["green"], label="Earth circles (negated)")
    ax.scatter(*arms[-1], color="k", s=8, zorder=5, label="chain tip (snapshot)")
    ax.scatter(0, 0, 0, color=OI["orange"], s=14, marker="*", zorder=5, label="Earth (origin)")
    lim = float(np.max(np.abs(np.vstack([truth, trail])))) * 1.02
    zl = max(float(np.max(np.abs(np.vstack([truth, trail])[:, 2]))) * 1.2, 0.3 * lim)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-zl, zl)
    ax.set_box_aspect((1, 1, zl / lim))  # true 1:1:1 scale
    ax.set_xlabel("x [AU]", labelpad=-2); ax.set_ylabel("y [AU]", labelpad=-2)
    ax.set_zlabel("z [AU]", labelpad=-2)
    ax.tick_params(pad=-2, labelsize=5.5)
    ax.zaxis.set_major_locator(plt.MaxNLocator(3))
    ax.view_init(elev=24, azim=-58)
    ax.legend(loc="upper left", fontsize=5.5, handlelength=1.2)
    save(fig, "fig_chain_3d.pdf")
    err = float(np.max(np.linalg.norm(trail - truth, axis=1)))
    return {"n_circles": len(chain), "synodic_days": float(syn), "tip_err_AU": err}


# ---------------------------------------------------------------- (e)
def load_study(path: Path = PERTURBED_JSON) -> dict:
    if not Path(path).exists():
        raise FileNotFoundError(f"{path} missing: run `uv run python experiments/perturbed_ceres.py` first")
    return json.loads(Path(path).read_text())


def fig_perturbed_ceres(d: dict) -> dict:
    s = d["study"]
    Ks = sorted(int(k) for k in s["naff_errors"])
    full = [s["naff_errors"][str(k)]["full_max"] for k in Ks]
    inter = [s["naff_errors"][str(k)]["interior_max"] for k in Ks]
    base = s["kepler_other"]
    fig, ax = plt.subplots(figsize=(W1, 2.7))
    ax.loglog(Ks, inter, "o-", ms=2.5, color=OI["blue"], label="NAFF, interior 80%")
    ax.loglog(Ks, full, "s-", ms=2.5, color=OI["sky"], label="NAFF, full window")
    ax.axhline(base["kepler_mean_a"]["full_max"], color=OI["verm"], ls="--", lw=0.9,
               label="Kepler, mean $a$")
    ax.axhline(base["kepler_fit6"]["interior_max"], color=OI["green"], ls="--", lw=0.9,
               label="Kepler, 6-element fit (interior)")
    ax.axhline(base["kepler_fit6"]["full_max"], color=OI["green"], ls=":", lw=0.9,
               label="Kepler, 6-element fit (full)")
    ax.set_xlabel(r"NAFF lines kept $K$ (circles $=2K+1$)")
    ax.set_ylabel("max geocentric error [AU] (200 yr)")
    ax.legend(loc="lower left", handlelength=1.6)
    save(fig, "fig_perturbed_ceres.pdf")
    osc = max(v["full_max"] for v in s["kepler_osculating_chain_errors"].values())
    return {"K": Ks, "full": full, "interior": inter, "mean_a": base["kepler_mean_a"]["full_max"],
            "fit6_int": base["kepler_fit6"]["interior_max"], "fit6_full": base["kepler_fit6"]["full_max"],
            "osc": osc, "span": s["span_years"]}


# ---------------------------------------------------------------- exports
def export_data(data_dir: Path = DATA_DIR, perturbed_json: Path = PERTURBED_JSON) -> list[dict]:
    """Write {mars,ceres,phaethon}_N{1,5,20,60}.json, perturbed_ceres.json and index.json."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    index = []
    for body in (MARS, CERES, PHAETHON):
        for N in (1, 5, 20, 60):
            chain = geocentric_chain(body, EARTH, N, N)
            meta = {"body": body.name, "earth": EARTH.name, "N_body": N, "N_earth": N,
                    "t_ref": body.epoch}
            name = f"{body.name.lower()}_N{N}.json"
            (data_dir / name).write_text(
                json.dumps(chain_to_json(chain, meta), separators=(",", ":")) + "\n")
            index.append({"file": name, "body": body.name, "N": N, "kind": "analytic",
                          "n_circles": len(chain)})
    d = load_study(perturbed_json)
    shutil.copyfile(perturbed_json, data_dir / "perturbed_ceres.json")
    index.append({"file": "perturbed_ceres.json", "body": d["meta"]["body"],
                  "N": int(d["study"]["chain_K"]), "kind": "naff",
                  "n_circles": len(d["epicycles"])})
    (data_dir / "index.json").write_text(json.dumps({"files": index}, indent=1) + "\n")
    return index


# ---------------------------------------------------------------- captions
def write_captions(a, b, c, dd, e):
    def f(x):
        return f"{x:.2g}"
    rows = ", ".join(f"$e={k}$: $\\rho={v:.3f}$" for k, v in a.items())
    capt = f"""# Figure captions (draft; numbers generated by experiments/make_figures.py)

## fig_coeff_decay.pdf
Magnitude of the Kepler-orbit Fourier coefficients $|c_k|$ (open markers, $k>0$, prograde) and
$|c_{{-k}}|$ (filled markers, retrograde; marker shape distinguishes $e$) in units of $a$, for $e\\in\\{{0.0934, 0.25, 0.6, 0.89\\}}$
(Mars, and three illustrative values; $e=0.89$ is Phaethon). Dashed lines: the Kapteyn rate
$\\rho(e)^k$, $\\rho=e\\,e^{{s}}/(1+s)$, $s=\\sqrt{{1-e^2}}$. {rows}. Decay is geometric
with ratio $\\rho$; retrograde terms are smaller than prograde ones by roughly $(e/(1+s))^2$ at small $k$.

## fig_geo_convergence.pdf
Maximum geocentric position error [AU] over a 50-year window (daily samples from JD 2461200.5)
of the truncated epicycle chain with $N_{{\\rm body}}=N_\\oplus=N$ versus exact Kepler motion,
for Mars, Ceres and Phaethon (solid), with the rigorous tail bound $\\epsilon_N({{\\rm body}})+\\epsilon_N(\\oplus)$
(dashed). Error first drops below $10^{{-6}}$ AU at $N={b['MARS']['N_1e6']}$ (Mars),
{b['CERES']['N_1e6']} (Ceres), {b['PHAETHON']['N_1e6']} (Phaethon), and below $10^{{-9}}$ AU at
{b['MARS']['N_1e9']}, {b['CERES']['N_1e9']}, {b['PHAETHON']['N_1e9']} (on the sampled grid of $N$).
Errors saturate at the double-precision roundoff floor ($\\sim2\\times10^{{-13}}$ AU).

## fig_mars_retrograde.pdf
Geocentric Mars over 2.75 years from 2020-05-31 (two oppositions). Top: sky track in ecliptic
longitude/latitude. Bottom: unwrapped longitude minus its mean motion ($0.524^\\circ$/day); grey bands
mark retrograde motion of the true path. Black: exact Kepler; dashed: Ptolemaic-order truncation
(Mars $k=+1$ deferent plus Earth $k=+1$ epicycle only, the model A of theory/derivation.md Part 4.3d); blue: $N=5$ (11+11
circles). The Ptolemaic model reproduces the retrograde loops but is displaced by up to
{c['window_ptol_max_deg']:.1f}$^\\circ$ in longitude in this window (planar 60-year sup norm
{c['planar_60yr']['max_deg']:.1f}$^\\circ$, RMS {c['planar_60yr']['rms_deg']:.1f}$^\\circ$, matching the
derivation's $33.2^\\circ$/$9.7^\\circ$); $N=5$ is within {c['window_N5_max_deg']:.3f}$^\\circ$.

## fig_chain_3d.pdf
Geocentric epicycle chain of Eros (inclination $10.8^\\circ$) with $N=3$ for Eros and for Earth
({dd['n_circles']} circles, each drawn in its own orbital plane about the preceding arm tip; red: body
circles, green: negated Earth circles). Blue: tip trail over one synodic period
({dd['synodic_days']:.0f} d); grey: exact path. Maximum tip error {dd['tip_err_AU']:.2f} AU at this low order.

## fig_perturbed_ceres.pdf
Maximum geocentric position error [AU] of a NAFF-extracted epicycle chain ($2K+1$ circles) fitted to
an N-body (Sun, Earth, Jupiter, Saturn) geocentric Ceres trajectory over {e['span']} years, versus the number of lines $K$:
interior 80% of the window (dark blue) and full window (light blue). Horizontal lines: best-fit
Kepler baselines: mean-$a$ Kepler {f(e['mean_a'])} AU; 6-element least-squares Kepler fit
{f(e['fit6_int'])} AU (interior), {f(e['fit6_full'])} AU (full). NAFF reaches {f(e['interior'][-1])} AU interior at
$K={e['K'][-1]}$ (full window {f(e['full'][-1])} AU; edge effects from the finite window grow with $K$),
beating the fitted Kepler orbit by a factor ~{e['fit6_int']/e['interior'][-1]:.0f} interior. The
osculating-element Kepler orbit (error {f(e['osc'])} AU) is omitted: that figure is an element-offset artifact,
not a fair comparison.
"""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    (FIG_DIR / "captions.md").write_text(capt)


def main():
    style()
    study = load_study()
    a = fig_coeff_decay()
    b = fig_geo_convergence()
    c = fig_mars_retrograde()
    dd = fig_chain_3d()
    e = fig_perturbed_ceres(study)
    write_captions(a, b, c, dd, e)
    idx = export_data()
    print(f"figures -> {FIG_DIR}; {len(idx)} data files -> {DATA_DIR}")


if __name__ == "__main__":
    main()
