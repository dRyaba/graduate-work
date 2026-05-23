# Cross-check v5 (2026-05-11) — pre-repetitions snapshot

First run with canonical m0..m5 labels and K4 excluded from the fleet,
but **before** the N=8 repetitions overhaul. Kept for traceability —
v6 is the citation-worthy run.

## Configuration

- 8 graphs: sausage 3×3 chains (3/4/5/6-block) + Geant2004 + Geant2009
  + IEEE-118 + UPS-Russia. K4 removed.
- 6 methods (m0..m5) with canonical labels.
- Per-(graph, method, d) **isolated** single calls — one
  `runWithTimeout` per row, no repetitions.
- Reverse method order within each cell (m5 → m4 → m3 → m2 → m1 → m0).
- Auto-ranged diameters: `d ∈ {dist, dist+1, dist+2, dist+3}`.
- Tiered timeouts: `kMethodTimeoutsSec = {60, 30, 120, 300, 300, 300}` s.

## Summary

- **Total runs**: 192 = 8 graphs × 4 d × 6 methods.
- **98 OK / 94 TIMEOUT / 0 ERROR**.
- **0 pair-disagreements > 1e-10**.
- Wall time: ~3 h (similar to v4).

## Why superseded

v6 (in `experiments/2026-05-12/`) adds:
- Three more sausage variants (`s3` slabs): 2-block, 4-block, 6-block.
- N=8 prepetitions per cell with median time + early-exit on two
  consecutive timeouts.
- `kMethodTimeoutsSec` raised to 1800 s for every method so m0/m1/m2
  also get a fair chance on the larger graphs.
- New CSV columns: `CompletedRuns`, `MinTimeSec`, `MaxTimeSec`.

v5 stays here only as the last single-run snapshot. Do not cite for
performance — its `TimeSec` values are single-trial measurements and
the m4/m5 sub-millisecond rows have high variance.

## Files

- `cross_check_v5.csv` — 192 raw single-call rows.
- `cross_check_v5.log` — stdout with per-cell ETA.
- `cross_check_agreement_matrix.csv` — 6×6 pair-agreement counts.
- `cross_check_disagreements.csv` — empty (no disagreements).
- `cross_check_trivial_cells.csv` — empty.
- `cross_check_perf_summary.csv` — per-method × per-family median/p95.
