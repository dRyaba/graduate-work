# Experiments — reproducing the cross-method consistency run

Goal: for every graph in the test grid, run all six reliability methods
(m0..m5) against the same `(s, t, d)` cell and check that every pair
agrees within `1e-10`. The output feeds the consistency chapter of the
thesis.

## Auto-ranging diameters

The CLI no longer accepts hand-picked diameters for `--cross-check`.
Instead, for each `(graph, s, t)` we compute `dist = dist_G(s, t)` via
Dijkstra and test

```
d ∈ { dist, dist + step, dist + 2·step, ..., dist + (count-1)·step }
```

clipped to `|V| - 1`. This avoids the degenerate regime `d < dist(s, t)`,
where every method trivially returns `R = 0` and "all methods agree" is a
tautology rather than an experiment.

Defaults: `--d-count 4 --d-step 1`. Tune via CLI.

## Running the full cross-check

```bash
cmake --build --preset release
./build/graph_reliability.exe --cross-check --output cross_check_results.csv 2>&1 \
    | tee cross_check.log
python scripts/cross_check_analyze.py cross_check_results.csv
```

The run takes ~25–35 minutes on an i7-class laptop under the tiered
timeout defaults. Progress prints per cell:

```
[3/25] 3_blocks_sausage_3x3_kao.txt d=11 s=0 t=28
    m0: TIMEOUT (>60s)
    m3: R=0.819793714074 (0.697s)
    ...
    [3/25 cells, elapsed=02:17, avg=45.7s/cell, ETA=16:40]
```

For a quick smoke test (K4 + sausage only, ~2 min):

```bash
./build/graph_reliability.exe --cross-check --quick --output quick.csv
```

## CSV column reference

`cross_check_results.csv` (one row per method × cell):

| Column | Meaning |
|---|---|
| `Graph` | Filename inside `graphs_data/` |
| `S`, `T`, `D` | Source vertex, target vertex, diameter bound |
| `MethodId` | 0..5 (see [ALGORITHMS.md](ALGORITHMS.md)) |
| `Method` | Human-readable method label |
| `Status` | `OK` / `TIMEOUT` / `ERROR` |
| `Reliability` | 15-digit `R(G, s, t, d)`; empty on non-OK |
| `TimeSec` | Wall time; empty on TIMEOUT |
| `Recursions` | Per-method counter (branching factor proxy) |
| `Error` | Exception message on ERROR |

## Tiered per-method timeouts

Default `--timeout 30` activates a tiered budget per method (in
`src/TestSuite.cpp`):

| Method | Timeout (s) | Rationale |
|---|---|---|
| m0 Standard Factoring | 60 | Baseline, should clear sausage-3 d=10 |
| m1 Recursive Decomposition | 30 | Known slow, don't spend more |
| m2 Simple Factoring | 120 | Mid-weight |
| m3 M-Decomposition | 300 | Reference workhorse |
| m4 Cancela-Petingi | 120 | Path-based, variance high |
| m5 M-Decomp + CPFM | 300 | Hybrid, matches m3 budget |

Passing any `--timeout N` other than `30` switches to uniform timeout `N`.

## Post-processing

`scripts/cross_check_analyze.py` (stdlib only):

```bash
python scripts/cross_check_analyze.py cross_check_results.csv
```

Outputs beside the input CSV:

- `cross_check_agreement_matrix.csv` — 6×6 counts `agreement[i][j]` =
  number of non-trivial cells where both `mᵢ` and `mⱼ` completed OK
  and `|Rᵢ − Rⱼ| ≤ tolerance`.
- `cross_check_disagreements.csv` — pair-disagreements above tolerance.
  Empty = methods agree.
- `cross_check_trivial_cells.csv` — cells where every OK method returned
  `R = 0` exactly. Should be empty after the auto-ranging fix; non-empty
  here signals a bug or a disconnected `(s, t)`.
- `cross_check_perf_summary.csv` — median and p95 `TimeSec` per method,
  globally and per graph family.

Flags:
- `--tolerance 1e-10` — override equality threshold.
- `--keep-trivial` — include `R = 0` cells in the agreement matrix
  (useful for diagnosing v1-style runs).

## Historical runs

- [experiments/2026-04-23/](../experiments/2026-04-23/) — v1 run with
  hardcoded diameters. Archived; **do not cite for agreement claims**
  (includes many trivial `R = 0` cells).

Later runs overwrite the repo-root CSV. Once a run is "final" for a
particular thesis chapter, move its artefacts to
`experiments/<YYYY-MM-DD>/` with a one-paragraph README describing what
changed since the previous archive.
