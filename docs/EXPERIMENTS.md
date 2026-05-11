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

## Per-(graph, method, d) isolated calls + reverse method order

Every row in the CSV is the result of **one independent call** to the
chosen method on the chosen `(s, t, d)` cell — `runCrossCheck` does not
share work between diameters even when the algorithm internally builds a
full `R(d)` curve. This matches the methodology of `run_benchmarks.sh` /
`benchmark_results.csv` and gives directly comparable per-d wall times
and recursion counts.

Within each `(graph, s, t, d)` cell, methods iterate in **reverse**
order m5 → m4 → m3 → m2 → m1 → m0. Rationale: m0/m1/m2 are the heavy
baselines that often time out on real-world graphs; running them last
means a `Ctrl+C` / TaskStop in the middle of a cell loses only the
heavy m0..m2 measurements while the m3/m4/m5 numbers that drive the
thesis chapter are already on disk.

## Tiered per-method timeouts

Default `--timeout 30` activates a tiered per-method budget defined in
`src/TestSuite.cpp` (`kMethodTimeoutsSec`):

| Method | Timeout (s) | Rationale |
|---|---|---|
| m0 Pure Factoring | 60 | Baseline; only K4 finishes within budget |
| m1 Block + Pure Facto | 30 | Diagnostic baseline, mostly TIMEOUT outside K4/sausage-3 |
| m2 Block + Conv + Simple Facto | 120 | Mid-weight; finishes most sausage cells |
| m3 Block + Conv + Modified Facto | 300 | Reference workhorse |
| m4 Cancela-Petingi (Nesterov) | 300 | Path-based, may enumerate many paths |
| m5 Block + Conv + Modified Cancela-Petingi | 300 | Pure m5 (no global fallback); same budget as m3 |

Passing any `--timeout N` other than `30` switches to a uniform timeout
`N` for every method (and ignores any `--method-timeout` overrides).

To patch the tiered defaults selectively, use `--method-timeout`:

```bash
./build/graph_reliability.exe --cross-check \
    --method-timeout m0=120,m1=60,m2=300 \
    --output cross_check.csv
```

Tokens are `mN=SEC` separated by commas. Unspecified methods keep their
tiered default. Useful e.g. when sampling sausage-3 carefully and you
want m0 to have a real chance.

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
