# Graph Reliability Analysis — Project Guide

C++17 library for **Diameter-Constrained s-t Reliability (DCR)** on
undirected graphs with independent edge failures. Implements six algorithms
(m0..m5) for computing `R(G, s, t, d) = Pr[∃ s-t path of length ≤ d]`
and cross-validates them against each other.

Context: graduate thesis of Novosibirsk State University (defense June 2026).

## Layout

- [include/graph_reliability/](include/graph_reliability/) — public headers.
  - `Graph.h` / `ReliabilityGraph.h` — core CSR graph + six algorithms.
  - `DataImporter.h` — KAO / edge-list parsers.
  - `TestSuite.h` — cross-check and benchmarking harness.
  - `CliHandlers.h` — dispatch table for subcommands.
- [src/](src/) — implementations mirroring each header.
  - [main.cpp](src/main.cpp) is now thin: global flag parse, dispatch.
  - [CliHandlers.cpp](src/CliHandlers.cpp) owns every `--subcommand`.
  - [ReliabilityGraph.cpp](src/ReliabilityGraph.cpp) (~1700 lines) is the
    algorithmic core — **do not refactor pre-defense**.
- [tests/](tests/) — Google Test (auto-fetched by CMake).
  - `test_ReliabilityGraph.cpp` includes the parametrised
    `MethodsAgreementTest` (regression guard for m0..m5 consistency).
- [graphs_data/](graphs_data/) — KAO/edge-list fixtures. See
  [graphs_data/README.md](graphs_data/README.md).
- [docs/](docs/) — thesis-oriented docs.
  - [docs/ALGORITHMS.md](docs/ALGORITHMS.md) — per-method theory and pseudocode.
  - [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md) — how to reproduce the cross-check.
  - [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), [docs/USER_GUIDE.md](docs/USER_GUIDE.md).
- [scripts/](scripts/) — Python post-processors (stdlib only, no pandas).
- [experiments/](experiments/) — archived run artefacts, one subdir per date.

## Build

MSYS2 MinGW64 toolchain is required (clean PATH: `export
PATH=/c/msys64/mingw64/bin:/usr/bin`). CMake presets in
[CMakePresets.json](CMakePresets.json):

```bash
cmake --preset mingw-release
cmake --build --preset release         # produces build/graph_reliability.exe
```

Run unit tests:

```bash
ctest --test-dir build --output-on-failure
```

## CLI quick start

```bash
./build/graph_reliability.exe --help
./build/graph_reliability.exe --run graphs_data/K4_kao.txt 0 3 3 4 1
./build/graph_reliability.exe --cross-check --output cross_check_results.csv
python scripts/cross_check_analyze.py cross_check_results.csv
```

## Conventions

- C++17, no exceptions across ABI boundaries. All errors via
  `graph_reliability::InvalidFormatException` family from `Exceptions.h`.
- Logger: `LOG_INFO(fmt, args...)` etc. from `Logger.h`. Macros no-op when
  `GRAPH_RELIABILITY_ENABLE_LOGGING` is off.
- Tests use Google Test with `INSTANTIATE_TEST_SUITE_P` for parametrisation.
- Algorithms mutate the graph in place; `TestSuite::runCrossCheck` keeps an
  in-memory template map and copies per run.

## Don't

- Don't refactor `ReliabilityGraph.cpp` before defense — cross-validation
  baseline changes are very expensive to re-audit.
- Don't commit files to the repo root that match
  `cross_check*.csv`, `*.log`, `*.svg`, `*.dot` — they're in `.gitignore`
  for a reason. Put them in `experiments/<date>/`.
- Don't widen the public API of `ReliabilityGraph` without updating
  `docs/ALGORITHMS.md` and the `MethodsAgreementTest` parameter list.
