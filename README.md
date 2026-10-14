# epicycles

Newtonian orbits as exact epicycle (Fourier) systems. (Full README: later bead.)

## Regenerate figures and SPA data

```bash
uv run python experiments/perturbed_ceres.py   # only if experiments/out/perturbed_ceres.json is missing (~3 min)
uv run python experiments/make_figures.py      # ~20 s: figures/*.pdf, figures/captions.md, web/public/data/*.json
uv run pytest -q
```
