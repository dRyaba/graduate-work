# Cross-check v7 (2026-05-22, Phase 2 of v6)

v7 is the citation-worthy follow-up to v6: same 11-graph fleet and
6-method matrix, but every cell that v6 finished with only one
single-shot timing was rerun with the full 8 repetitions so that
**TimeSec is a true median**, not a single-trial point.

## What changed vs v6

- Started from v6 (`experiments/2026-05-22/cross_check_v6.csv`),
  patched the CSV in place: dropped the 49 `Status=OK,CompletedRuns=1`
  rows (single-shot screening cells), kept everything else (157
  full-8-rep rows + 56 TIMEOUTs + 2 ERRORs).
- New column `TimeoutBudgetSec` recorded per-row (m0/m1 = 300 s,
  m2..m5 = 1800 s). Disambiguates TIMEOUT rows — previously the
  budget was implicit.
- Rerun (`--cross-check --resume --repetitions 8 --method-timeout
  m0=300,m1=300`) added 49 fresh rows back in, each with full 8-rep
  median, min, max.
- Adaptive cap was active (≥30 min first OK → 1 rep; 5..30 → 2;
  <5 → 8), so the budget held.

## Phase 2 summary

```
Summary: total=49 OK=49 TIMEOUT=0 ERROR=0
No discrepancies > 0.0
All methods agree within tolerance
```

Wall time of Phase 2 itself: **1 h 39 min**.

## v7 overall (post-merge: v6 carry-over + Phase 2)

- **Total rows**: 264 = 11 graphs × 4 d × 6 methods.
- **0 pair-disagreements > 1e-10**.
- **0 trivial cells**.
- **2 ERROR** rows (UPS-Russia d=18, d=19, m3 `std::bad_alloc`).
- **CompletedRuns=8**: 168 rows.
- **CompletedRuns=4**: 1 row (adaptive cap).
- **CompletedRuns=0** (TIMEOUT/ERROR): 95 rows.
- **CompletedRuns=1**: **0 rows** — single-shot screening fully replaced.

## Agreement matrix (excl. trivial)

```
      m0  m1  m2  m3  m4  m5
m0     4   4   4   4   4   4
m1     4  24  24  24  19  24
m2     4  24  32  32  22  32
m3     4  24  32  35  25  35
m4     4  19  22  25  33  31
m5     4  24  32  35  31  41
```

m5 finishes 41/44 cells (only 3 TIMEOUTs on UPS-Russia). m3 35 + 2
ERR. m4 33. m2 32. m1 24. m0 4 (sausage-2 s3 only).

## Per-family completion

| Family | m0 | m1 | m2 | m3 | m4 | m5 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| sausage 3×3 (16) | 0/16 | 5/16 | 16/16 | 16/16 | 9/16 | 16/16 |
| sausage s3 (12) | 4/12 | 16 | 12/12 | 12/12 | 9/12 | 12/12 |
| Geant2004 (4) | 0/4 | 3/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| Geant2009 (4) | 0/4 | 0/4 | 0/4 | 0/4 | 4/4 | 4/4 |
| IEEE-118 (4) | 0/4 | 0/4 | 0/4 | 3/4 | 4/4 | 4/4 |
| UPS-Russia (4) | 0/4 | 0/4 | 0/4 | 0+2 ERR | 3/4 | 1/4 |

m4 + m5 are the only methods that finish every Geant2009 / IEEE-118
cell. UPS-Russia is the hardest of the fleet.

## Per-family perf (median / p95 sec, OK rows only)

| Family | m1 | m2 | m3 | m4 | m5 |
|---|---:|---:|---:|---:|---:|
| sausage | 8.34 / 121.8 | 0.64 / 9.85 | 0.90 / 13.49 | 0.47 / 439.3 | 0.17 / 8.55 |
| Geant2004 | 94.85 / 158.8 | 41.40 / 67.99 | 0.010 / 0.16 | 0.0002 / 0.001 | 0.0008 / 0.001 |
| Geant2009 | — | — | — | 0.0004 / 0.002 | 0.0065 / 0.018 |
| IEEE-118 | — | — | 37.89 / 700.0 | 0.0003 / 0.006 | 0.0037 / 0.436 |
| UPS-Russia | — | — | — | 0.14 / 4.82 | 0.001 / 0.001 |

## Files

- `cross_check_v7.csv` — 264 rows, header now 15 columns incl.
  `TimeoutBudgetSec`.
- `cross_check_v7.log` — phase 2 stdout (resume + 49 new cells).
- `cross_check_agreement_matrix.csv`, `cross_check_disagreements.csv`,
  `cross_check_trivial_cells.csv`, `cross_check_perf_summary.csv` —
  analyzer outputs.

## Reproducing

```bash
# Phase 0 + Phase 1 (v6, ~13 h, archived in experiments/2026-05-22/)
./build/graph_reliability.exe --cross-check \
    --method-timeout m0=300,m1=300 \
    --output cross_check_v6.csv

# Patch the v6 CSV in place: drop CompletedRuns=1 OK rows and add
# the TimeoutBudgetSec column (the awk one-liner is in
# experiments/2026-05-22-v7/README.md commit message).

# Phase 2 (~1.5 h)
./build/graph_reliability.exe --cross-check --resume \
    --method-timeout m0=300,m1=300 \
    --repetitions 8 \
    --output cross_check_v7.csv

python scripts/cross_check_analyze.py cross_check_v7.csv
```

This is the run to cite in the thesis comparison chapter.
