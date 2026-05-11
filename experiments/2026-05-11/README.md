# Cross-check v4 (2026-05-11) — full 6-method per-d isolated run

First clean run after the methodology overhaul:
- Every CSV row is **one isolated call** to `(graph, method, d)` —
  no shared work between diameters, matches `run_benchmarks.sh` /
  `benchmark_results.csv`. Each row's `TimeSec` and `Recursions` are
  per-call, not amortized.
- Methods iterate in **reverse** within each cell (m5 → m4 → m3 →
  m2 → m1 → m0). A `Ctrl+C` mid-cell would lose only heavy m0/m1/m2
  rows; the m3/m4/m5 numbers that drive the thesis chapter are
  written first.
- m5 has **no fallback** to m4 anymore (removed 2026-05-02). When a
  block exceeds the per-block CPFM threshold, m5 honestly times out
  instead of silently delegating to global CPFM.

## Configuration

```
./build/graph_reliability.exe --cross-check --output cross_check_v4.csv
```

- `d_count=4 d_step=1` (defaults), so each (graph, s, t) tests
  `d ∈ {dist, dist+1, dist+2, dist+3}` clipped to |V|-1.
- Tiered timeouts: m0=60, m1=30, m2=120, m3=300, m4=300, m5=300 s.
- Methods: 0, 1, 2, 3, 4, 5 (full set).

## Summary

- **Total runs**: 210 = 9 graphs × 4 d × 6 methods − 6 (K4 only has
  3 d's: dist=1 → {1, 2, 3}, clipped at |V|-1=3).
- **115 OK / 95 TIMEOUT / 0 ERROR**.
- **0 pair-disagreements > 1e-10** — every pair of methods that both
  finished OK on the same cell agreed within tolerance.
- Wall time: **3:14:07**.

## Per-method completion

| Method | OK | TIMEOUT | Notes |
|---|---:|---:|---|
| **m5** M-Decomp + CPFM | 32/35 | 3 | Only TIMEOUTs on UPS-Russia |
| **m4** Cancela-Petingi | 27/35 | 8 | TIMEOUTs on dense sausage cells (d ≥ 16) and one UPS-Russia |
| **m3** M-Decomposition | 25/35 | 10 | TIMEOUTs on Geant2009 + half of IEEE-118 + UPS-Russia |
| **m2** Simple Factoring | 21/35 | 14 | Same pattern as m3 |
| **m1** Recursive Decomposition | 7/35 | 28 | Only K4 + part of sausage-3 |
| **m0** Standard Factoring | 3/35 | 32 | Only K4 — fundamentally exponential without decomposition |

Note: m0 TIMEOUT on sausage matches v2 baseline; not a regression.
m0 was not exercised on sausage in `run_benchmarks.sh`
(`run_benchmarks.sh` loops `m=3,4,5` only) so the legacy CSV did
not have m0 sausage data either.

## Per-family breakdown

| Family | m0 | m1 | m2 | m3 | m4 | m5 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| K4 (3 d) | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 |
| sausage 3×3 (16 cells: 4 graphs × 4 d) | 0/16 | 4/16 | 16/16 | 16/16 | 9/16 | 16/16 |
| Geant2004 (4 d) | 0/4 | 0/4 | 2/4 | 4/4 | 4/4 | 4/4 |
| Geant2009 (4 d) | 0/4 | 0/4 | 0/4 | 0/4 | 4/4 | 4/4 |
| IEEE-118 (4 d) | 0/4 | 0/4 | 0/4 | 2/4 | 4/4 | 4/4 |
| UPS-Russia (4 d) | 0/4 | 0/4 | 0/4 | 0/4 | 3/4 | 1/4 |

m4 and m5 are the only methods that finish every Geant2009 and
IEEE-118 cell within budget. UPS-Russia is the hardest: m5 manages
only 1/4 (the smallest d), m4 manages 3/4 — without the old fallback,
m5 honestly hits the per-block CPFM wall on this graph's biconnected
core.

## Files

- `cross_check_v4.csv` — 211 lines (210 data + header).
- `cross_check_v4.log` — full stdout incl. per-cell ETA.
- `cross_check_agreement_matrix.csv` — 6×6 pair-agreement counts
  (excluding trivial cells, but there are 0 trivial cells anyway).
- `cross_check_disagreements.csv` — empty (no disagreements).
- `cross_check_trivial_cells.csv` — empty.
- `cross_check_perf_summary.csv` — median / p95 per method, ALL +
  per family (K4, sausage, Geant2004, Geant2009, IEEE-118, UPS-Russia).

## Reproducing

```bash
./build/graph_reliability.exe --cross-check --output cross_check_v4.csv \
    2>&1 | tee cross_check_v4.log
python scripts/cross_check_analyze.py cross_check_v4.csv
```

Or via the launcher:

```bash
RUN_DIR=experiments/$(date +%F) METHODS=0,1,2,3,4,5 \
    scripts/run_cross_check.sh
```
