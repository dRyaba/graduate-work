# User Guide

## Building the Project

### Prerequisites

- **C++17** compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake 3.16+**
- **OpenMP** (optional, for parallel processing)
- **Doxygen** (optional, for documentation)

### Build Steps

#### Debug Build (with logging)

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build .
```

#### Release Build (optimized, no logging)

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

#### Debug Build with Sanitizers

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON ..
cmake --build .
```

## Running the Application

### Command Line Interface

#### Run Single Test

```bash
./graph_reliability --run [--1-based] <file> <s> <t> <d> <method> [reps]
```

**Parameters:**

- `--1-based`: (optional) If present, s and t are 1-based (as in VKR/course work documents)
- `file`: Graph file path (KAO format)
- `s`: Source vertex (0-based, or -1 to auto-detect)
- `t`: Target vertex (0-based, or -1 to auto-detect)
- `d`: Diameter upper bound
- `method`: Method ID (0-4)
- `reps`: Number of repetitions (optional, default: 1)

**Example:**

```bash
./graph_reliability --run graphs_data/K4_kao.txt 0 3 2 3
./graph_reliability --run --1-based graphs_data/Geant2004_kao.txt 12 96 8 3  # 1-based (VKR)
./graph_reliability --run graphs_data/3_blocks_sausage_3x3_kao.txt -1 -1 10 3 5
```

#### Run Comprehensive Tests

```bash
./graph_reliability --test <method_id> [output_file]
```

**Example:**

```bash
./graph_reliability --test 3 results.csv
```

#### Convert Formats

**Edge List to KAO:**

```bash
./graph_reliability --convert edge2kao input.edgelist output.kao 0.9
```

**KAO to Edge List:**

```bash
./graph_reliability --convert kao2edge input.kao output.edgelist
```

#### Verbose Mode (Debug builds only)

```bash
./graph_reliability --verbose --run graphs_data/K4_kao.txt 0 3 2 3
```

## Method Selection

| ID | Method | Description | Best For |
|----|--------|-------------|----------|
| 0 | Standard Factoring | Baseline, no decomposition | Small graphs |
| 1 | Recursive Decomposition | Nested recursion | Academic comparison |
| 2 | Simple Factoring | Convolution + simple factoring | Medium graphs |
| 3 | M-Decomposition | Convolution + modified factoring | Large graphs (recommended) |
| 4 | Cancela-Petingi | Path-based factoring with SPT | Alternative, decomposed structures |

## Input Data Formats

### KAO Format

```
0,2,5,8,10,13,17,21,24,27,31,35,38,40,43,46,48
2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,13,6,9,11,14,7,10,12,15,8,11,16,9,14,10,13,15,11,14,16,12,15
0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
0.9
```

**Format:**

- Line 1: KAO offsets (comma-separated)
- Line 2: FO adjacency list (comma-separated, 1-based)
- Line 3: Target vertices (comma-separated, with leading 0)
- Line 4: Edge reliability (single value or per-edge)

### Edge List Format

```
1 -- 2
2 -- 3
3 -- 4
// Comments are supported
```

**Format:**

- One edge per line: `u -- v`
- Comments start with `//`
- Vertices are 1-based

## Output Interpretation

### Single Test Output

```
Running test on: graphs_data/K4_kao.txt
  Source: 0, Target: 3
  Diameter: 2, Method: M-Decomposition (Level 3)
  Repetitions: 1

Results:
  Reliability: 0.810000
  Avg Time:    0.001234 seconds
  Avg Recs:    15
```

**Fields:**

- **Reliability**: Probability of path existence (0.0 to 1.0)
- **Avg Time**: Average execution time in seconds
- **Avg Recs**: Average number of recursive calls

### CSV Test Results

```csv
Graph,S,T,D,Reps,Method,AvgTime,AvgRel,AvgRecs
K4_kao.txt,0,3,2,1,M-Decomposition (Level 3),0.001234,0.810000,15
```

## Running Tests

### Unit Tests

```bash
cd build
cmake --build . --target graph_reliability_tests
./tests/graph_reliability_tests
```

Or use CMake test target:

```bash
cd build
ctest
```

### Integration Tests

```bash
cd build
cmake --build . --target run_tests
```

## Troubleshooting

### Common Issues

#### "File not found" error

- Ensure `graphs_data/` directory exists
- Check file paths are correct
- Use absolute paths if relative paths fail

#### "Invalid format" error

- Verify KAO format is correct (4 lines)
- Check Edge List format (one edge per line)
- Ensure no extra whitespace

#### Slow performance

- Use Method 3 (M-Decomposition) for large graphs
- Reduce diameter constraint if possible
- Use Release build for better performance
- Check if graph can be decomposed into smaller blocks

#### High memory usage

- Large graphs may require significant memory
- Consider using block decomposition
- Reduce number of repetitions

### Debug Mode

Enable verbose logging in Debug builds:

```bash
./graph_reliability --verbose --run <file> <s> <t> <d> <method>
```

Check log file: `graph_reliability.log`

### Performance Analysis

1. Compare methods: Run same graph with different methods
2. Check recursion count: Higher count = more computation
3. Monitor execution time: Use Release build for accurate timing
4. Profile with tools: Use `perf`, `valgrind`, or `gprof`

## Best Practices

1. **Use Method 3** for production calculations
2. **Release builds** for performance-critical runs
3. **Debug builds** for development and debugging
4. **Multiple repetitions** for statistical accuracy
5. **Validate input** before running large calculations
6. **Monitor memory** for very large graphs

## Examples

### Example 1: Small Graph

```bash
# Calculate reliability for K4 graph
./graph_reliability --run graphs_data/K4_kao.txt 0 3 2 3
```

### Example 2: Large Graph with Multiple Repetitions

```bash
# Run 5 repetitions for statistical accuracy
./graph_reliability --run graphs_data/6_blocks_sausage_3x3_kao.txt -1 -1 24 3 5
```

### Example 3: Comprehensive Testing

```bash
# Run all test graphs with M-Decomposition
./graph_reliability --test 3 comprehensive_results.csv
```

### Example 4: Format Conversion

```bash
# Convert Edge List to KAO format
./graph_reliability --convert edge2kao my_graph.edgelist my_graph.kao 0.9

# Use converted file
./graph_reliability --run my_graph.kao 0 5 10 3
```
