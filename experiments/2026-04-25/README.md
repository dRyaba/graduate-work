# Cross-check v3 (2026-04-25) — m3/m4/m5 with grouped multi-d API

First targeted m3-vs-m4-vs-m5 run after the multi-diameter CDF APIs
landed. Each (graph, s, t) group is now factored once per method via
`calculateReliabilityCdf*` and one row is emitted per d ∈ [d_min,
d_max]. Time per row is the **per-(s, t) wall time**, not a fictitious
per-d slice.

## Summary

- **Cells**: 27 non-trivial (K4 × 3d, sausage 3×3 chains × 4d each,
  Geant2004 × 4d, IEEE-118 × 4d).
- **Runs**: 27 × 3 = 81. **61 OK / 20 TIMEOUT / 0 ERROR**.
- **Pair-disagreements > 1e-10**: **0**.
- **Trivial cells**: 0.
- Invocation: `--cross-check --methods 3,4,5` (default tiered timeouts;
  m4 bumped to 300 s to match m3/m5 since multi-d CPFM enumerates all
  paths up to d_max in one call).

## Agreement matrix (`cross_check_agreement_matrix.csv`)

```
      m3  m4  m5
m3    23   7  23
m4     7  11  11
m5    23  11  27
```

m5 clears every cell (27/27); m3 25/27 (single-cell IEEE-118 timeouts);
m4 11/27 — multi-d CPFM is dominated by path enumeration up to d_max
and times out on every sausage chain at d_max ≥ 13.

## Per-family performance (median / p95 seconds, OK rows)

- **K4 (3d)**: all three methods sub-millisecond.
- **sausage (3×3 chains, 16 cells)**: m3 median 1.93 s / p95 3.22 s,
  m5 1.52 s / 2.50 s, m4 entirely TIMEOUT.
- **Geant2004 (4d)**: m3 0.40 s, m5 3.8 ms, m4 21.5 ms — m5 fastest,
  m4 close behind, m3 noticeably slower.
- **IEEE-118 (4d)**: m3 TIMEOUT × 4, m4 46.0 s, m5 30.3 s.

## Reading the time column

Every row in `cross_check_v3.csv` for methods 3/4/5 carries the same
`TimeSec` value across all d in its (graph, s, t) group: that value is
the wall time of **one** factorization that produced the entire d-curve.
Per-d times are not separable for these methods. The `Recursions`
column is non-zero only on the first row of each group and reports the
total across all d.

## Reproducing

```bash
./build/graph_reliability.exe --cross-check --methods 3,4,5 \
    --output cross_check_v3.csv 2>&1 | tee cross_check_v3.log
python scripts/cross_check_analyze.py experiments/2026-04-25/cross_check_v3.csv
```
