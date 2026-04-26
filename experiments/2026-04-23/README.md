# Cross-check v1 artifacts (2026-04-23)

First full run of `--cross-check` before the auto-ranging fix. Kept as a
historical snapshot; **do not cite for agreement claims**.

## Why archived

The diameter list was hardcoded per graph (e.g. `{9, 10, 11, 12}` for
`3_blocks_sausage_3x3_kao.txt`). Several of those diameters are below
`dist(s, t)`, producing trivially `R = 0` for every method. Counting those
rows as "6/6 methods agree within 1e-10" is tautological, not a real
cross-check.

After this run the CLI was changed to auto-range diameters from `dist(s, t)`
upward (`--d-count N --d-step N`), so every cell is guaranteed
`R ∈ (0, R_classical]`. See v2 results in the repo root.

## Contents

- `cross_check_results.csv` — 150 raw (graph, s, t, d, method) rows.
- `cross_check.log` — stdout with ETA progress.
- `cross_check_agreement_matrix.csv` — 6×6 agreement counts (inflated by
  trivial cells).
- `cross_check_disagreements.csv` — empty.
- `cross_check_perf_summary.csv` — median / p95 timings per method × family.
- `*.svg`, `*.dot`, `graph_widget.html` — visualizer outputs generated
  during this period.

## Reproducing the analysis (not the run)

```bash
python scripts/cross_check_analyze.py --keep-trivial experiments/2026-04-23/cross_check_results.csv
```

`--keep-trivial` preserves the v1 matrix shape; without it the Python
post-processor strips the all-zero cells and reports only
`len(non_trivial)` cells.
