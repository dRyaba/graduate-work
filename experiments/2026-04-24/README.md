# Cross-check v2 (2026-04-24) — auto-ranged diameters

First run after the methodological fix: diameters for every `(graph, s, t)`
cell are auto-ranged from `dist(s, t)` upward
(`--d-count 4 --d-step 1`), so no cell is in the degenerate `d < dist`
regime. Every OK row has `R > 0`.

## Summary

- **Cells**: 27 non-trivial (K4 × 3d, sausage 3×3-{3,4,5,6}block × 4d,
  sausage-3 4×4 absent here, Geant2004 × 4d, IEEE-118 × 4d).
- **Runs**: 27 × 6 = 162. **102 OK / 60 TIMEOUT / 0 ERROR**.
- **Pair-disagreements > 1e-10**: **0** — all methods agree on every
  cell where both complete OK.
- **Trivial cells (all methods return R = 0)**: **0** — methodology fix
  confirmed.
- **Wall time**: 84:50.

## Agreement matrix (`cross_check_agreement_matrix.csv`)

```
      m0  m1  m2  m3  m4  m5
m0     3   3   3   3   3   3
m1     3   8   8   8   8   8
m2     3   8  20  20  12  20
m3     3   8  20  25  17  25
m4     3   8  12  17  19  19
m5     3   8  20  25  19  27
```

Diagonal = cells where method completed OK. m5 clears every cell (27/27),
m3 clears 25/27, m2 20/27, m4 19/27, m1 8/27, m0 3/27.

## Per-family performance (median / p95 seconds, OK rows only)

- **K4**: all methods under 1 ms.
- **sausage (3×3 chains)**: m5 fastest (median 1.2 s), m3 close behind
  (median 2.5 s), m4 highly variable (p95 67 s on saturated cases), m1
  bottlenecks (only 5/16 OK).
- **Geant2004**: m2 median 53 s (one cell OK), m3/m4/m5 sub-second or
  single-ms. m0/m1 never complete.
- **IEEE-118**: m3 completes 2/4 (p95 252 s), m4/m5 always under 12 ms.
  m0/m1/m2 entirely TIMEOUT.

## Contents

- `cross_check_v2.csv` — 162 raw rows.
- `cross_check_v2.log` — stdout including the cell plan, per-cell ETAs.
- `cross_check_agreement_matrix.csv` — the 6×6 table above.
- `cross_check_disagreements.csv` — empty (no disagreements > 1e-10).
- `cross_check_trivial_cells.csv` — empty (no R = 0 cells).
- `cross_check_perf_summary.csv` — median / p95 per method × family.

## Reproducing the analysis

```bash
python scripts/cross_check_analyze.py experiments/2026-04-24/cross_check_v2.csv
```
