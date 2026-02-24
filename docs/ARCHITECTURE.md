# Architecture Documentation

## Overview

Graph Reliability Analysis Library is a C++17 library for calculating network reliability with diameter constraints. The library provides efficient algorithms for reliability calculations using various decomposition and factoring methods.

## Project Structure

```
graduate-work/
├── include/                    # Header files
│   ├── graph_reliability.h    # Main header
│   └── graph_reliability/     # Library headers
│       ├── Graph.h            # Base graph structure (CSR format)
│       ├── ReliabilityGraph.h # Extended graph with reliability methods
│       ├── DataImporter.h     # Data loading and conversion
│       ├── TestSuite.h        # Testing framework
│       ├── GraphOperations.h  # Graph utilities
│       ├── Logger.h           # Logging system
│       └── Exceptions.h       # Custom exceptions
├── src/                       # Source files
│   ├── main.cpp              # CLI application
│   ├── Graph.cpp             # Graph implementation
│   ├── ReliabilityGraph.cpp  # Reliability calculations
│   ├── DataImporter.cpp      # Data import implementation
│   ├── TestSuite.cpp         # Test suite implementation
│   ├── GraphOperations.cpp   # Graph utilities
│   └── Logger.cpp            # Logging implementation
├── tests/                     # Unit tests
│   ├── CMakeLists.txt        # Test build configuration
│   ├── test_main.cpp         # Test entry point
│   ├── test_Graph.cpp        # Graph tests
│   ├── test_ReliabilityGraph.cpp # Reliability tests
│   ├── test_DataImporter.cpp # DataImporter tests
│   ├── test_GraphOperations.cpp # GraphOperations tests
│   └── test_config.json      # Test configuration
├── graphs_data/              # Test data files
├── docs/                     # Documentation
├── CMakeLists.txt           # Build configuration
└── README.md               # Project readme
```

## Core Components

### Graph Class

Base graph structure using **Compressed Sparse Row (CSR)** format for efficient memory usage:

- **KAO (Offset Array)**: Points to the start of each vertex's adjacency list in FO
- **FO (Flat Array)**: Adjacency list containing all neighbors
- **PArray**: Probability array for each edge

**Key Operations:**
- `numVertices()`, `numEdges()`: Graph statistics
- `hasEdge()`, `findEdge()`: Edge queries
- `calculateDistance()`: Shortest path using Dijkstra
- `removeEdge()`, `swapVertices()`: Graph modifications
- `findArticulationPoints()`, `decomposeByVertices()`: Graph analysis

### ReliabilityGraph Class

Extends `Graph` with target vertices and reliability calculation methods:

**Reliability Calculation Methods:**
1. **Standard Factoring** (Level 0): Baseline method, no decomposition
2. **Recursive Decomposition** (Level 1): Nested recursion (inefficient)
3. **Simple Factoring** (Level 2): Convolution with simple factoring
4. **M-Decomposition** (Level 3): Convolution with modified factoring (fastest)

**Key Features:**
- Target vertices for reliability calculations
- Multiple diameter constraint algorithms
- Block decomposition for large graphs
- Parallel processing support (OpenMP)

### DataImporter Class

Centralized data management:

- Loads graphs from KAO and Edge List formats
- Converts between formats
- Manages data file paths
- Validates loaded graphs

### Logger System

Conditional compilation logging:

- **Debug builds**: Full logging with spdlog (file + console)
- **Release builds**: Zero overhead (all logging removed by compiler)
- Periodic progress logging for long-running calculations
- Configurable log levels

### Exception Hierarchy

Custom exceptions for better error handling:

```
GraphReliabilityException (base)
├── InvalidGraphException
├── FileNotFoundException
├── InvalidFormatException
├── InvalidParameterException
└── CalculationException
```

## Data Formats

### KAO Format

```
0,2,5,8,10,13,17,21,24,27,31,35,38,40,43,46,48  # KAO offsets
2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,...          # FO adjacency list
0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0                # Target vertices
0.9                                              # Edge reliability
```

### Edge List Format

```
1 -- 2
2 -- 3
3 -- 4
// Comments are supported
```

## Build System

### CMake Configuration

- **Debug**: Logging enabled, debug symbols, no optimization
- **Release**: No logging, full optimization, zero overhead
- **RelWithDebInfo**: Optional logging, optimized with debug info
- **MinSizeRel**: Size-optimized build

### Dependencies

- **spdlog**: Logging (only in Debug builds)
- **Google Test**: Unit testing framework
- **OpenMP**: Parallel processing (optional)
- **C++17**: Required standard

## Performance Considerations

1. **CSR Format**: Efficient memory usage for sparse graphs
2. **Conditional Logging**: Zero overhead in Release builds
3. **Block Decomposition**: Reduces complexity for large graphs
4. **M-Decomposition**: Fastest algorithm for reliability calculations
5. **OpenMP**: Parallel processing for CPU-intensive operations

## Testing

- **Unit Tests**: Google Test framework for individual components
- **Integration Tests**: CLI-based tests for end-to-end scenarios
- **Test Configuration**: JSON-based test configuration

## Logging Strategy

- **Key Points Only**: Start/end of calculations, not every recursion
- **Periodic Progress**: Every 100k recursions (configurable)
- **Error Logging**: All errors and warnings
- **Release Builds**: Complete removal of logging code
