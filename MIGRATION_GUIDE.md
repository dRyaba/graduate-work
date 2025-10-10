# Migration Guide: From Legacy to Professional Architecture

## Overview

This guide explains how to migrate from the legacy codebase to the new professional architecture while preserving all core functionality.

## Key Changes

### 1. File Structure

**Before:**
```
graduate-work/
├── graduate work.cpp          # ❌ Space in filename
├── kGraphOperations.cpp       # ❌ Mixed naming
├── GraphOperations.cpp        # ❌ Global namespace
└── DataImporter.cpp           # ❌ No organization
```

**After:**
```
graduate-work/
├── include/                   # ✅ Organized headers
│   ├── graph_reliability.h   # ✅ Main header
│   └── graph_reliability/    # ✅ Namespace organization
├── src/                      # ✅ Source organization
│   └── main.cpp             # ✅ Clean main file
├── CMakeLists.txt           # ✅ Professional build system
└── README.md                # ✅ Professional documentation
```

### 2. Naming Conventions

**Before:**
```cpp
struct kGraph { ... };           // ❌ Inconsistent naming
void Factoring(...);             // ❌ Global functions
int Nconst;                      // ❌ Global variables
std::vector<int> cutPoints;      // ❌ Global state
```

**After:**
```cpp
namespace graph_reliability {
    class ReliabilityGraph { ... };  // ✅ Consistent naming
    void performFactoring(...);      // ✅ Class methods
    // ✅ No global variables
}
```

### 3. Class Design

**Before:**
```cpp
// ❌ Mixed responsibilities
struct Graph {
    std::vector<int> KAO, FO;
    std::vector<double> PArray;
    // Mixed with reliability calculations
};

struct kGraph : public Graph {
    std::vector<int> Targets;
    // Reliability-specific methods mixed in
};
```

**After:**
```cpp
// ✅ Clear separation of concerns
class Graph {
    // Pure graph structure and operations
};

class ReliabilityGraph : public Graph {
    // Reliability-specific functionality
};
```

### 4. Error Handling

**Before:**
```cpp
// ❌ Inconsistent error handling
std::ifstream fin("path");
if (!fin) {
    std::cout << "Error!\n";
    throw std::runtime_error("OPEN_ERROR");
}
```

**After:**
```cpp
// ✅ Consistent exception handling
try {
    auto graph = importer.loadKAOGraph(filename);
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
}
```

### 5. Memory Management

**Before:**
```cpp
// ❌ Manual memory management
kGraph G = kGraphFileInput(fin);
// Potential memory leaks
```

**After:**
```cpp
// ✅ Smart pointers
auto graph = importer.loadKAOGraph(filename);
// Automatic memory management
```

## Migration Steps

### Step 1: Update Includes

**Before:**
```cpp
#include "kGraphOperations.h"
#include "DataImporter.h"
```

**After:**
```cpp
#include "graph_reliability.h"
using namespace graph_reliability;
```

### Step 2: Update Class Names

**Before:**
```cpp
kGraph G = kGraphFileInput(fin);
G.ReliabilityDiamConstr2VertMDecompose(x, y, d);
```

**After:**
```cpp
auto graph = importer.loadKAOGraph(filename);
auto result = graph->calculateReliabilityWithMDecomposition(x, y, d);
```

### Step 3: Update Method Calls

**Before:**
```cpp
// Global functions
Factoring(G, variant, d, reliability);
GraphMerging(k);
```

**After:**
```cpp
// Class methods
graph->performFactoring(variant, d, reliability);
GraphOperations::mergeGraphsFromFile(k, importer);
```

### Step 4: Update Data Loading

**Before:**
```cpp
std::ifstream fin("C:/Users/User/source/repos/graduate work/graduate work/input.txt");
kGraph G = kGraphFileInput(fin);
```

**After:**
```cpp
DataImporter importer("graphs_data/");
auto graph = importer.loadKAOGraph("input.txt");
```

## Backward Compatibility

The new architecture maintains full backward compatibility for:
- All reliability calculation algorithms
- Graph data formats (KAO, Edge List)
- Test configurations and results
- Performance characteristics

## Benefits of Migration

1. **Professional Structure**: Clean, organized codebase
2. **Better Maintainability**: Clear separation of concerns
3. **Improved Testing**: Comprehensive test framework
4. **Cross-Platform**: Works on all major platforms
5. **Modern C++**: Uses C++17 features and best practices
6. **Documentation**: Professional API documentation
7. **Build System**: CMake-based professional build system

## Testing Migration

To verify the migration:

```bash
# Build new version
mkdir build && cd build
cmake .. && make

# Run tests
./graph_reliability --test 0

# Compare results with legacy version
diff test_results.csv ../test_summary_results.csv
```

## Support

For migration questions or issues:
- Check the API documentation
- Review the example code in `src/main.cpp`
- Open an issue for specific problems
