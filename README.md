# Graph Reliability Analysis Library

![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)

A comprehensive C++ library for network reliability analysis, providing efficient algorithms for calculating network reliability with diameter constraints and various decomposition methods.

## Features

- **Multiple Reliability Algorithms**: Support for M-decomposition, recursive decomposition, simple factoring, and standard factoring
- **Flexible Graph Representation**: CSR (Compressed Sparse Row) format for memory efficiency
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Professional API**: Clean, well-documented C++17 interface
- **Comprehensive Testing**: Unit tests (Google Test) and integration tests with configurable parameters
- **Format Conversion**: Convert between Edge List and KAO formats
- **Parallel Processing**: OpenMP support for performance-critical operations
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
1. Установите MSYS2: https://www.msys2.org/
2. В MSYS2 выполните: `pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake`
3. Добавьте `C:\msys64\mingw64\bin` в PATH

Подробные инструкции см. в **[SETUP.md](SETUP.md)**

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

```bash
# Convert Edge List to KAO format
./graph_reliability --convert edge2kao input.edgelist output.kao 0.9

# Convert KAO to Edge List format
./graph_reliability --convert kao2edge input.kao output.edgelist

# Run comprehensive tests
./graph_reliability --test 0 results.csv

# Show help
./graph_reliability --help
```

## Project Structure

```plaintext
graduate-work/
├── include/                    # Header files
│   ├── graph_reliability.h    # Main header
│   └── graph_reliability/     # Library headers
│       ├── Graph.h            # Base graph structure
│       ├── ReliabilityGraph.h # Extended graph for reliability
│       ├── DataImporter.h     # Data import system
│       ├── TestSuite.h        # Testing framework
│       └── GraphOperations.h  # Graph utilities
├── src/                       # Source files
│   ├── main.cpp              # Main application
│   ├── Graph.cpp             # Base graph implementation
│   ├── ReliabilityGraph.cpp  # Reliability calculations
│   ├── DataImporter.cpp      # Data import implementation
│   ├── TestSuite.cpp         # Testing framework
│   └── GraphOperations.cpp   # Graph utilities
├── graphs_data/              # Test data (KAO and EdgeList files)
├── docs/                     # Documentation (Doxygen)
├── build/                    # Build output (gitignored)
├── CMakeLists.txt           # Build configuration
└── README.md                # This file
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

Comprehensive testing framework for different reliability calculation methods.

### Reliability Calculation Methods

1. **M-Decomposition** (`method_id = 0`): Advanced decomposition with M-factorization
2. **Recursive Decomposition** (`method_id = 1`): Recursive decomposition with simple factoring
3. **Simple Factoring** (`method_id = 2`): Basic factoring with decomposition
4. **Standard Factoring** (`method_id = 3`): Traditional factoring method

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

The library is optimized for performance with:

- CSR format for memory efficiency
- OpenMP parallelization for CPU-intensive operations
- Efficient algorithms with complexity analysis
- Minimal memory allocations

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

## Documentation

### User Documentation

- **[User Guide](docs/USER_GUIDE.md)**: Complete guide for using the library
- **[Architecture](docs/ARCHITECTURE.md)**: System architecture and design
- **[Algorithms](docs/ALGORITHMS.md)**: Detailed algorithm descriptions

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
- **OpenMP** (optional, for parallel processing)
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
