# Cross-check v6 (started 2026-05-12, finished 2026-05-22)

First end-to-end run with the full methodology overhaul described in
Part 4 of the plan. Spans 10 days because of two restarts (s3 marker
fix; m0/m1 timeout reduction; phased reduction in repetitions).

## Configuration

- 11 graphs:
  - sausage 3×3 chains: 3/4/5/6-block,
  - sausage s3 slabs:   2/4/6-block,
  - real-world:         Geant2004 + Geant2009 + IEEE-118 + UPS-Russia.
  - K4 excluded (too trivial, no signal).
- 6 methods (m0..m5) with canonical labels.
- Per-(graph, method, d) **isolated** single-d calls.
- Reverse method order m5 → m4 → m3 → m2 → m1 → m0 within each cell.
- Auto-ranged diameters: `d ∈ {dist, dist+1, dist+2, dist+3}`.
- Tiered timeouts: m0 = m1 = **300 s** (5 min, dropped from the
  default 1800 s after the heavy m0 cells started eating an hour
  per timeout); m2 = m3 = m4 = m5 = **1800 s** (30 min).
- **N=8 repetitions** with median time per cell, early-exit on
  two consecutive TIMEOUTs, adaptive cap (≥30 min first OK → 1 rep;
  5..30 min → 2 reps; <5 min → full 8).
- 157 cells from the early phase already had 8 reps in the CSV when
  the run was paused; the rest of the 264 total cells were finished
  with `--repetitions 1` in a "phase 1 screening" pass for time.

## Summary

- **Total runs**: 264 = 11 graphs × 4 d × 6 methods.
- 49 OK + 56 TIMEOUT + 2 ERROR in the final phase-1 pass; combined
  with the resumed 157 OK rows from the early phase, the CSV has
  data for every cell.
- **0 pair-disagreements > 1e-10** — methods agree across the entire
  fleet wherever they both finish.
- **0 trivial cells** (no R = 0 cells; auto-ranging worked).
- Wall time: 16:21 for the final phase-1 pass on top of ~14 h of
  the earlier resumed prefix.
- **ERROR rows**: 2 — UPS-Russia d=18 and d=19, method m3, both
  `std::bad_alloc` (modified factoring exhausting memory on the
  graph's dense biconnected core). Not a correctness issue;
  other methods on the same cells produce R, agreement is intact.

## Agreement matrix (excl. trivial cells, 44 non-trivial)

```
      m0  m1  m2  m3  m4  m5
m0     4   4   4   4   4   4
m1     4  24  24  24  19  24
m2     4  24  32  32  22  32
m3     4  24  32  35  25  35
m4     4  19  22  25  33  31
m5     4  24  32  35  31  41
```

m5 finishes 41/44 cells (only 3 TIMEOUTs on UPS-Russia). m3 manages
35 + 2 ERR ≈ 37; m4 33; m2 32; m1 24; m0 4 (sausage-2 s3 only).

## Per-family completion

| Family | Cells | m0 | m1 | m2 | m3 | m4 | m5 |
|---|---:|:-:|:-:|:-:|:-:|:-:|:-:|
| sausage 3×3 (16) | 16 | 0 | 5 | 16 | 16 | 9 | 16 |
| sausage s3 (12) | 12 | 4 | 16 | 12 | 12 | 9 | 12 |
| Geant2004 (4) | 4 | 0 | 3 | 4 | 4 | 4 | 4 |
| Geant2009 (4) | 4 | 0 | 0 | 0 | 0 | 4 | 4 |
| IEEE-118 (4) | 4 | 0 | 0 | 0 | 3 | 4 | 4 |
| UPS-Russia (4) | 4 | 0 | 0 | 0 | 0 (+2 ERR) | 3 | 1 |

m4 + m5 are the only methods that finish every Geant2009 / IEEE-118
cell within budget. UPS-Russia is the hardest of the fleet.

## Caveat on TimeSec

Cells from the early phase (≈ first 157 rows of the CSV) carry
medians over 1..8 repetitions per the adaptive policy. Cells from
the phase-1 screening pass carry single-trial wall times
(`CompletedRuns = 1`, `MinTimeSec == TimeSec == MaxTimeSec`). When
plotting / citing TimeSec, filter by `CompletedRuns` to avoid mixing
the two regimes:

```bash
awk -F, 'NR==1 || $11>=2' cross_check_v6.csv > v6_multi_run.csv
awk -F, 'NR==1 || $11==1' cross_check_v6.csv > v6_single_run.csv
```

A follow-up phase 2 (rerun all cells with `CompletedRuns == 1 &&
TimeSec < 300` at `--repetitions 8`) is a natural next step but is
not part of this run.

## Files

- `cross_check_v6.csv` — 264 rows.
- `cross_check_v6.log` — full stdout incl. ETA + per-cell breakdown.
- `cross_check_agreement_matrix.csv` — 6×6 pair-agreement counts.
- `cross_check_disagreements.csv` — empty (no disagreements).
- `cross_check_trivial_cells.csv` — empty.
- `cross_check_perf_summary.csv` — median / p95 per method × family.

## Reproducing

```bash
./build/graph_reliability.exe --cross-check \
    --method-timeout m0=300,m1=300 \
    --output cross_check_v6.csv 2>&1 | tee cross_check_v6.log
python scripts/cross_check_analyze.py cross_check_v6.csv
```

Or via the launcher: `scripts/run_cross_check.sh`.
