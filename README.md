# Graph Reliability Analysis Library

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/your-repo/graph-reliability)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

A comprehensive C++ library for network reliability analysis, providing efficient algorithms for calculating network reliability with diameter constraints and various decomposition methods.

## Features

- **Multiple Reliability Algorithms**: Support for M-decomposition, recursive decomposition, simple factoring, and standard factoring
- **Flexible Graph Representation**: CSR (Compressed Sparse Row) format for memory efficiency
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Professional API**: Clean, well-documented C++17 interface
- **Comprehensive Testing**: Built-in test suite with configurable parameters
- **Format Conversion**: Convert between Edge List and KAO formats
- **Parallel Processing**: OpenMP support for performance-critical operations

## Quick Start

### Building the Project

```bash
# Clone the repository
git clone <repository-url>
cd graduate-work

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the project
make -j$(nproc)

# Run tests
make run_tests
```

### Basic Usage

```cpp
#include "graph_reliability.h"

using namespace graph_reliability;

// Create data importer
DataImporter importer("graphs_data/");

// Load a graph
auto graph = importer.loadKAOGraph("example.kao");

// Calculate reliability
auto result = graph->calculateReliabilityWithMDecomposition(1, 10, 15);

std::cout << "Reliability: " << result.reliability << std::endl;
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

```
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
│   └── main.cpp              # Main application
├── graphs_data/              # Test data
├── docs/                     # Documentation
├── CMakeLists.txt           # Build configuration
└── README.md               # This file
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
```
0,2,5,8,10,13,17,21,24,27,31,35,38,40,43,46,48
2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,13,6,9,11,14,7,10,12,15,8,11,16,9,14,10,13,15,11,14,16,12,15
0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
0.9
// Comments are supported
```

### Edge List Format
```
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

Run the comprehensive test suite:

```bash
# Run all tests with method 0 (M-Decomposition)
./graph_reliability --test 0

# Run tests with custom output file
./graph_reliability --test 1 my_results.csv
```

## Documentation

Generate API documentation:

```bash
# Build documentation (requires Doxygen)
make docs

# View documentation
open docs/html/index.html
```

## Requirements

- **C++17** compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake 3.16+**
- **OpenMP** (optional, for parallel processing)
- **Doxygen** (optional, for documentation)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Based on research in network reliability analysis
- Implements advanced graph decomposition algorithms
- Optimized for real-world network analysis scenarios

## Support

For questions, issues, or contributions, please:
- Open an issue on GitHub
- Contact the development team
- Check the documentation in the `docs/` directory
