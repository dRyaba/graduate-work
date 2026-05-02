# Graph Reliability Analysis Library

![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)

A comprehensive C++ library for network reliability analysis, providing efficient algorithms for calculating network reliability with diameter constraints and various decomposition methods.

## Features

- **Multiple Reliability Algorithms**: Support for M-decomposition, recursive decomposition, simple factoring, standard factoring, and Cancela-Petingi path-based method
- **Flexible Graph Representation**: CSR (Compressed Sparse Row) format for memory efficiency
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Professional API**: Clean, well-documented C++17 interface
- **Comprehensive Testing**: Unit tests (Google Test) and integration tests with configurable parameters
- **Format Conversion**: Convert between Edge List and KAO formats
- **Logging System**: Conditional compilation logging (zero overhead in Release builds)
- **Custom Exceptions**: Hierarchical exception system for better error handling
- **Extensive Documentation**: Architecture, algorithms, and user guides

## Quick Start

### Building the Project

#### Предварительная проверка окружения

**Windows:**

```powershell
# Проверьте наличие всех инструментов
.\check_environment.ps1
```

**Linux/Mac:**

```bash
# Проверьте наличие компилятора и CMake
g++ --version
cmake --version
```

#### Установка зависимостей (если нужно)

**Windows (рекомендуется MinGW-w64):**

1. Установите MSYS2: <https://www.msys2.org/>
2. В MSYS2 выполните: `pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake`
3. Добавьте `C:\msys64\mingw64\bin` в PATH

Подробные инструкции см. в **[SETUP.md](SETUP.md)**

#### Recommended: CMake presets

The project ships [CMakePresets.json](CMakePresets.json) with `mingw-release`
and `mingw-debug` configurations pinned to the MSYS2 toolchain. Single
command from the repo root:

```bash
cmake --preset mingw-release
cmake --build --preset release         # produces build/graph_reliability.exe
ctest --test-dir build --output-on-failure
```

For a Debug build with logging enabled, swap `mingw-release` → `mingw-debug`
and `--preset release` → `--preset debug`.

#### Manual configuration (alternative)

If presets are not available (older CMake, non-MSYS2 toolchain), the manual
flow below still works.

#### Debug Build (with logging)

**Windows (MinGW):**

```powershell
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "MinGW Makefiles"
cmake --build . -j8
```

**Linux/Mac:**

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . -j$(nproc)
```

#### Release Build (optimized, no logging overhead)

**Windows (MinGW):**

```powershell
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
cmake --build . -j8
```

**Linux/Mac:**

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)
```

#### Debug Build with Sanitizers

**Linux/Mac:**

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON ..
cmake --build . -j$(nproc)
```

**Windows:** Санитайзеры работают только с Clang на Windows

### Basic Usage

```cpp
#include "graph_reliability.h"
#include <iomanip>
#include <iostream>

using namespace graph_reliability;

// Create data importer
DataImporter importer("graphs_data/");

// Load a graph
auto graph = importer.loadKAOGraph("K4_kao.txt");

// Calculate reliability
auto result = graph->calculateReliabilityWithMDecomposition(1, 10, 15);

std::cout << "Reliability: " << std::fixed << std::setprecision(15) << result.reliability << std::endl;
std::cout << "Execution time: " << result.execution_time_sec << " seconds" << std::endl;
```

### Command Line Interface

```powershell
# From project root (Windows)
.\build\graph_reliability.exe --help

# Run reliability calculation
.\build\graph_reliability.exe --run graphs_data\K4_kao.txt 0 3 3 4 1

# Convert Edge List to KAO format
.\build\graph_reliability.exe --convert edge2kao input.edgelist output.kao 0.9

# Convert KAO to Edge List format
.\build\graph_reliability.exe --convert kao2edge input.kao output.edgelist

# Run comprehensive tests
.\build\graph_reliability.exe --test 3 results.csv

# Cross-check all 6 methods against each other (diameters auto-ranged from dist(s,t))
.\build\graph_reliability.exe --cross-check --output cross_check_results.csv
python scripts\cross_check_analyze.py cross_check_results.csv
```

See [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md) for full cross-check options
(`--d-count`, `--d-step`, `--quick`, `--timeout`, `--keep-trivial` for the
post-processor) and how to reproduce the thesis consistency tables.

On Linux/macOS:

```bash
./build/graph_reliability --help
```

## Project Structure

```plaintext
graduate-work/
├── include/
│   ├── graph_reliability.h            # Umbrella header
│   └── graph_reliability/
│       ├── Graph.h                    # Base CSR graph
│       ├── ReliabilityGraph.h         # m0..m5 + multi-d CDF APIs
│       ├── DataImporter.h             # KAO / edge-list parsers
│       ├── GraphOperations.h          # Format conversion utilities
│       ├── TestSuite.h                # Cross-check + bench harness
│       ├── CliHandlers.h              # Subcommand dispatch
│       ├── GraphVisualizer.h          # FR layout, SVG/DOT export
│       ├── PathEnumerator.h           # Bounded-length s-t path enum
│       ├── CancelaPetingiState.h      # CPFM solver state
│       ├── SeriesParallelTransform.h  # Series/parallel reductions
│       ├── Logger.h                   # spdlog wrapper macros
│       └── Exceptions.h               # InvalidFormatException family
├── src/                               # Implementations mirroring headers
│   ├── main.cpp                       # Entry point + dispatcher
│   └── CliHandlers.cpp                # All --subcommand handlers
├── tests/                             # Google Test (auto-fetched)
│   ├── test_ReliabilityGraph.cpp      # incl. MethodsAgreementTest
│   └── test_config.json
├── graphs_data/                       # KAO / edge-list fixtures
│   └── README.md                      # Catalogue with sources
├── docs/
│   ├── ALGORITHMS.md                  # Per-method theory
│   ├── ARCHITECTURE.md
│   ├── EXPERIMENTS.md                 # Cross-check reproduction
│   └── USER_GUIDE.md
├── scripts/
│   └── cross_check_analyze.py         # stdlib post-processor
├── experiments/                       # Archived run artefacts (per date)
├── build/                             # gitignored
├── CMakeLists.txt
├── CMakePresets.json                  # mingw-release / mingw-debug
├── CLAUDE.md                          # Project navigation map
├── SETUP.md
└── README.md                          # This file
```

## API Reference

### Core Classes

#### `Graph`

Base graph structure using CSR format for efficient memory usage and fast adjacency queries.

#### `ReliabilityGraph`

Extended graph structure with target vertices and specialized reliability calculation methods.

#### `DataImporter`

Centralized system for loading graphs from various formats and managing data file paths.

#### `TestSuite`

Comprehensive testing framework for different reliability calculation
methods. Hosts `runCrossCheck` (auto-ranged diameters, per-cell timeout,
optional method filter) and the static helpers `chooseDiameters`,
`runWithTimeout`, `runWithTimeoutCdf`.

#### `GraphOperations`

Stateless utilities for KAO ↔ edge-list conversion and other graph
transformations not tied to a specific instance.

#### `GraphVisualizer`

Two-level Fruchterman-Reingold layout + SVG/DOT export. Powers the
`--visualize` subcommand.

#### `PathEnumerator`

Enumerates simple s-t paths up to a length bound; used by m4/m5.

#### `CancelaPetingiState`

State machine for path-based factoring (CPFM): ESS / ISPT / GlobalISPT
inclusion-exclusion bookkeeping.

#### `SeriesParallelTransform`

Series and parallel edge reductions used by the decomposition pipeline.

### Reliability Calculation Methods

| Method ID | Name | Description | Speed |
|-----------|------|-------------|-------|
| 0 | **Standard Factoring** | Baseline: direct edge factoring without decomposition. R = p×R(contract) + (1-p)×R(delete) | Slowest |
| 1 | **Recursive Decomposition** | Block decomposition with nested recursion. Academic comparison only | Very slow |
| 2 | **Simple Factoring** | Block decomposition + convolution + simple factoring | Fast |
| 3 | **M-Decomposition** | Block decomposition + modified factoring (computes all diameters in one pass) | Fastest |
| 4 | **Cancela-Petingi** | Path-based factoring with SPT/ISPT optimizations. Operates on path lists, not graphs | Efficient |
| 5 | **M-Decomp + CPFM** | Hybrid: block decomposition + multi-diameter path-based factoring inside blocks | Efficient |

**Recommended**: Use Method 3 (M-Decomposition) for production, Method 5 (M-Decomp + CPFM) for graphs with many paths where ISPT is effective.

#### Multi-diameter (CDF) APIs

For workflows that need `R(s, t, d)` for many diameters at once
(parameter sweeps, cross-check, plotting), m3 / m4 / m5 also expose a
single-pass CDF entry point that factors once over `[dist(s, t), d_max]`
and returns the entire vector:

```cpp
ReliabilityCdfResult cdf = graph->calculateReliabilityCdfMDecompositionCPFM(s, t, d_max);
// cdf.cdf[d] == R(G, s, t, d) for d ∈ [0, d_max]; zero for d < dist(s, t).
// cdf.execution_time_sec is the wall time of one factorization.
```

Available methods: `calculateReliabilityCdfMDecomposition` (m3),
`calculateReliabilityCdfCancelaPetingi` (m4),
`calculateReliabilityCdfMDecompositionCPFM` (m5). The single-d entry
points keep their existing signature.

See [docs/ALGORITHMS.md](docs/ALGORITHMS.md) for detailed algorithm descriptions.

## Data Formats

### KAO Format

```plaintext
0,2,5,8,10,13,17,21,24,27,31,35,38,40,43,46,48
2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,13,6,9,11,14,7,10,12,15,8,11,16,9,14,10,13,15,11,14,16,12,15
0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
0.9
// Comments are supported
```

### Edge List Format

```plaintext
1 -- 2
2 -- 3
3 -- 4
// Comments are supported
```

## Performance

### Benchmark Results (April 2026, Release build, -O3)

#### Sausage chain graphs at d = d_min

| Graph | d | m3 | m4 | **m5** | Speedup m5/m3 |
|---|---|---|---|---|---|
| 2-block 3×3 | 8 | 16 ms | 37 ms | 2.6 ms | ×6 |
| 3-block 3×3 | 10 | 581 ms | 0.12 ms | **0.20 ms** | **×2905** |
| 3-block 4×4 | 13 | TIMEOUT | 0.86 ms | **0.86 ms** | — |
| 4-block 3×3 | 14 | 1084 ms | 0.91 ms | **0.24 ms** | **×4517** |
| 5-block 3×3 | 18 | 1626 ms | 14.7 ms | **0.26 ms** | **×6254** |
| 6-block 3×3 | 22 | 2161 ms | 220 ms | **0.34 ms** | **×6356** |

#### Real-world networks at d = d_min

| Graph | V | E | d | m3 | **m4** | m5 |
|---|---|---|---|---|---|---|
| GEANT 2004 | 103 | 127 | 6 | 3 ms | **0.07 ms** | 0.70 ms |
| GEANT 2009 | 390 | 503 | 12 | TIMEOUT | **0.07 ms** | 0.06 ms |
| IEEE 118-node | 118 | 168 | 10 | 36.8 s | **0.08 ms** | 0.06 ms |
| UPS Russia | 63 | 108 | 18 | TIMEOUT | 158 ms | 160 ms |

All methods produce identical reliability values within `|R_i − R_j| ≤ 1e-10`,
empirically verified on 27 non-trivial (graph, s, t, d) cells across
K4, sausage chains, Geant2004 and IEEE-118 — see
[experiments/2026-04-25/](experiments/2026-04-25/README.md) for the
agreement matrix and per-(s, t) wall times under the multi-d CDF API.

**Key takeaway**: m5 (M-Decomp + CPFM) is optimal for structured sausage chains (up to ×6356 vs m3). m4 (Cancela-Petingi) is optimal for real-world networks with large biconnected components.

See [docs/ALGORITHMS.md](docs/ALGORITHMS.md#experimental-results) for full
analysis and [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md) for how to
reproduce the cross-check tables.

### Implementation optimizations

- CSR format for memory efficiency
- Gap optimization in m5: each block only computes CDF for its relevant diameter range
- ISPT / ESS / GlobalISPT path pruning in CPFM (m4, m5)
- Multi-diameter CDF entry points (m3 / m4 / m5) factor once per (s, t) and return the full R(d) curve

## Testing

### Unit Tests

Run unit tests using Google Test:

```bash
cd build
cmake --build . --target graph_reliability_tests
./tests/graph_reliability_tests

# Or use CTest
ctest
```

### Integration Tests

Run the comprehensive test suite:

```bash
# Run all tests with method 3 (M-Decomposition)
./graph_reliability --test 3 results.csv

# Run tests with method 0 (Standard Factoring)
./graph_reliability --test 0

# Run tests with custom output file
./graph_reliability --test 2 my_results.csv
```

### Test Configuration

Test configurations can be customized in `tests/test_config.json`.

### Cross-method consistency regression

`tests/test_ReliabilityGraph.cpp` includes a parametrised
`MethodsAgreementTest` that runs m0..m5 on the same (graph, s, t, d)
cells (K4 and sausage-3 cases inside CI-friendly per-method timeouts)
and asserts pair-wise agreement within `1e-10`. It is registered with
ctest as `MethodsAgreementTests` and runs together with the rest.

For a full cross-method consistency run with the CDF-based grouped
multi-d API:

```bash
./build/graph_reliability.exe --cross-check --output cross_check.csv \
    2>&1 | tee cross_check.log
python scripts/cross_check_analyze.py cross_check.csv
```

See [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md) for the methodology
(auto-ranged diameters, tiered timeouts, `--methods` filter) and the
post-processor's outputs.

## Documentation

### User Documentation

- **[CLAUDE.md](CLAUDE.md)**: Project navigation map (layout, conventions, don'ts)
- **[User Guide](docs/USER_GUIDE.md)**: Complete guide for using the library
- **[Architecture](docs/ARCHITECTURE.md)**: System architecture and design
- **[Algorithms](docs/ALGORITHMS.md)**: Detailed algorithm descriptions
- **[Experiments](docs/EXPERIMENTS.md)**: Reproducing the cross-method consistency run
- **[Graph catalogue](graphs_data/README.md)**: Source and size of every test graph

### API Documentation

Generate API documentation:

```bash
# Build documentation (requires Doxygen)
cmake --build . --target docs

# View documentation
open docs/html/index.html
```

### Debugging Tips

- Use Debug builds with `--verbose` flag for detailed logging
- Check `graph_reliability.log` for log output (Debug builds only)
- Use sanitizers for memory error detection: `-DENABLE_SANITIZERS=ON`
- Release builds have zero logging overhead for production use

## Requirements

- **C++17** compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake 3.16+**
- **Doxygen** (optional, for documentation)
- **spdlog** (automatically fetched via CMake, only for Debug builds)
- **Google Test** (automatically fetched via CMake, for unit tests)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add unit tests for new functionality
5. Ensure all tests pass (`ctest`)
6. Update documentation if needed
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

### Code Style

- Follow C++17 best practices
- Use const correctness
- Prefer RAII and smart pointers
- Add logging for key operations (Debug builds only)
- Use custom exceptions for error handling
- Write unit tests for new features

## Acknowledgments

- Based on research in network reliability analysis
- Implements advanced graph decomposition algorithms
- Optimized for real-world network analysis scenarios

## Support

For questions, issues, or contributions, please:

- Open an issue on GitHub
- Contact the development team
- Check the documentation in the `docs/` directory
