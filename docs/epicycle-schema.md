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
| `source` | `"body"` or `"earth"` | Which orbit the circle came from. Earth circles are already negated (phase shifted by pi); do not negate again. |
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
