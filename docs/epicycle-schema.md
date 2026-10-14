# Epicycle chain JSON schema (version 1)

Produced by `epicycles.export.write_chain` / `chain_to_json`; consumed by the SPA.
Test `tests/test_expansion.py::test_export_round_trip` checks that evaluating this
JSON with the formula below reproduces `evaluate_chain`.

```json
{
  "schema_version": 1,
  "meta": {"body": "Mars", "earth": "Earth", "N_body": 12, "N_earth": 6,
           "units": {"length": "AU", "time": "day", "angle": "rad"},
           "t_ref": 2451595.0},
  "epicycles": [
    {"radius": 1.52, "omega": 0.0091, "phase": 0.3,
     "normal": [0,0,1], "u": [1,0,0], "v": [0,1,0], "source": "body", "k": 1}
  ]
}
```

## Fields

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | int | Always `1`. |
| `meta.body`, `meta.earth` | string / null | Names of the observed body and of the observer's planet. |
| `meta.N_body`, `meta.N_earth` | int / null | Truncation orders: harmonics `|k| <= N` kept for each. |
| `meta.units` | object | `length: "AU"`, `time: "day"`, `angle: "rad"` (fixed). |
| `meta.t_ref` | number | Julian date (TDB, days) at which every `phase` is given. |
| `epicycles` | array | Circles, sorted by `radius` descending (outermost/largest first). |
| `radius` | number >= 0 | Circle radius [AU]. |
| `omega` | number | Signed angular rate [rad/day]. Positive = counter-clockwise about `+normal`; negative = clockwise; `0` = fixed arm (the k=0 offset). |
| `phase` | number | Angle [rad] from `u` toward `v` at `t = t_ref`. Not wrapped to a fixed interval. |
| `normal` | [x,y,z] | Unit normal of the circle's plane, ecliptic J2000 frame. |
| `u`, `v` | [x,y,z] | Orthonormal in-plane basis, `u x v = normal`. |
| `source` | `"body"`, `"earth"` or `"naff"` | Which orbit the circle came from. Earth circles are already negated (phase shifted by pi); do not negate again. `"naff"` marks circles extracted numerically by NAFF frequency analysis from a (perturbed) trajectory; these are already geocentric, so `meta.earth`/`N_*` may be null and `k` is null (loaders must not assume `body`/`earth` only). |
| `k` | int / null | Harmonic index: circle comes from `c_k exp(i k M)` (k may be negative = retrograde). `null` if not from a Kepler expansion. |

## Evaluation formula (exact)

The geocentric ecliptic J2000 position at Julian date `t` is the plain sum over all circles:

```
theta_j(t) = phase_j + omega_j * (t - t_ref)
r(t)       = sum_j radius_j * ( cos(theta_j) * u_j + sin(theta_j) * v_j )      [AU]
```

For drawing nested arms in list order, the j-th arm tip is the partial sum of the
first j terms (origin = Earth); the final tip equals `r(t)`.

## Notes

* All circles are centred on the origin (Earth); there is no per-circle centre offset.
  Nesting is just the order in which terms are added.
* Truncation error is bounded by `epicycles.bounds.geocentric_tail_bound`
  (theory/derivation.md, Part 4); the JSON stores only the truncated terms.
* Body and Earth orbital planes differ, so `normal` generally varies between circles.

## Extension: `study` key (perturbed_ceres.json)

`web/public/data/perturbed_ceres.json` (copied from `experiments/out/perturbed_ceres.json`) is a normal
schema-v1 chain (all circles `source = "naff"`, 2K+1 = 401 circles for K = 200) plus one optional
top-level key `study` (ignored by chain evaluators; added keys are backwards compatible):

| `study` field | Meaning |
|---|---|
| `masses_Msun`, `epoch_jd`, `span_years`, `dt_days`, `n_samples` | N-body setup and sampling. |
| `window`, `interior` | NAFF window description; "interior" = middle 80% of samples. |
| `K_max`, `chain_K` | Lines extracted / lines in the stored chain (200). |
| `naff_errors` | `{K: {n_circles, full_max, full_rms, interior_max, interior_rms}}` [AU] vs the N-body truncated to the first K lines. |
| `kepler_osculating_chain_errors` | Same error dict for an osculating-element Kepler chain (dominated by element offset; not a fair comparison). |
| `kepler_other` | Fitted Kepler baselines `kepler_mean_a`, `kepler_fit6` (with `a`, `e`). |
| `kepler_note`, `units` | Free text; `units` = `AU`. |

`web/public/data/index.json`: `{"files": [{"file", "body", "N", "kind": "analytic"|"naff", "n_circles"}, ...]}`
(for the naff entry, `N` is `chain_K`).
