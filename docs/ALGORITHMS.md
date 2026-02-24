# Algorithms Documentation

## Overview

This document describes the reliability calculation algorithms implemented in the Graph Reliability Analysis Library.

## Problem Definition

Given a graph G(V, E) with edge reliabilities p(e) ∈ [0, 1], source vertex s, target vertex t, and diameter constraint d, calculate the probability that there exists a path from s to t with length ≤ d.

## Algorithm Methods

### Method 0: Standard Factoring (Level 0)

**Description**: Baseline method using direct factoring without decomposition.

**Algorithm**:
1. Select an unreliable edge e
2. Factor: R = p(e) × R(e=reliable) + (1-p(e)) × R(e=removed)
3. Recursively solve subproblems

**Complexity**: Exponential in number of edges

**Use Case**: Baseline comparison, small graphs

### Method 1: Recursive Decomposition (Level 1)

**Description**: Uses true nested recursion with block decomposition.

**Algorithm**:
1. Decompose graph into k-blocks
2. For each block i:
   - For each path length L in block i:
     - Recursively call solver for block i+1
3. Combine results

**Complexity**: Very high due to nested recursion

**Use Case**: Academic comparison, demonstrates inefficiency

### Method 2: Simple Factoring with Decomposition (Level 2)

**Description**: Uses block decomposition with convolution and simple factoring.

**Algorithm**:
1. Decompose graph into k-blocks
2. Calculate minimum diameter for each block
3. Compute reliability for each block for diameter range [min, min+gap]
4. Convolve block reliabilities using Migov's formula
5. Use simple factoring within blocks

**Complexity**: O(B × E × D) where B=blocks, E=edges, D=diameter range

**Use Case**: Medium-sized graphs, good balance

### Method 3: M-Decomposition (Level 3)

**Description**: Fastest method using modified factoring with block decomposition.

**Algorithm**:
1. Decompose graph into k-blocks
2. Build block graph
3. Find path through block graph
4. For each block in path:
   - Calculate reliability for ALL diameters in one pass (modified factoring)
   - Returns cumulative distribution function (CDF)
5. Convolve block reliabilities efficiently

**Key Optimization**: Modified factoring computes reliabilities for multiple diameters simultaneously, avoiding redundant calculations.

**Complexity**: O(B × E × D) with better constant factors

**Use Case**: Large graphs, production use

## Block Decomposition

### K-Blocks

A k-block is a maximal 2-connected subgraph. The graph is decomposed into k-blocks connected by articulation points.

**Algorithm**:
1. Find articulation points using DFS
2. Identify 2-connected components
3. Build block graph (blocks as nodes, articulation points as edges)

### Block Graph Path Finding

Find shortest path through block graph from source block to target block.

## Factoring Methods

### Standard Factoring

For edge e with probability p:
- R = p × R(e=reliable) + (1-p) × R(e=removed)

### Modified Factoring

Computes reliability for multiple diameters simultaneously:
- Maintains CDF: P(distance ≤ d) for d ∈ [lower, upper]
- More efficient than computing each diameter separately

## Distance Calculation

### Dijkstra's Algorithm

Used for shortest path calculations:
- Priority queue implementation
- Early termination when distance exceeds constraint
- Complexity: O(E log V)

### Floyd-Warshall

Used for distance constraint checking:
- All-pairs shortest paths
- Complexity: O(V³)

## Performance Characteristics

| Method | Complexity | Best For |
|--------|-----------|----------|
| Standard Factoring | O(2^E) | Small graphs (< 20 edges) |
| Recursive Decomposition | O(B × 2^E) | Academic comparison |
| Simple Factoring | O(B × E × D) | Medium graphs |
| M-Decomposition | O(B × E × D) | Large graphs |

## Implementation Details

### Recursion Counting

All methods track recursion count for performance analysis.

### Time Measurement

Execution time measured using high-resolution clock.

### Memory Management

- CSR format minimizes memory usage
- Smart pointers for automatic memory management
- Move semantics to avoid unnecessary copies

## Optimization Techniques

1. **Early Termination**: Stop when distance constraint violated
2. **Memoization**: Cache results for repeated subproblems (future enhancement)
3. **Block Decomposition**: Reduce problem size
4. **Modified Factoring**: Compute multiple diameters simultaneously
5. **Parallel Processing**: OpenMP for independent calculations

## Scientific Sources

The algorithms are based on the following works:

- **Рябинин Д.В.** ВКР бакалавр. Разработка метода расчёта надёжности сети с ограничением на диаметр. (Algorithm 1: Modified Factoring, Algorithm 2: ReliabilityDiamConstr2VertMDecompose)
- **Рябинин Д.В.** Курсовая работа (магистратура). Декомпозиционный подход к расчёту надёжности сети с ограничением на диаметр. (Обобщённая формула (2.3), итеративная свёртка)
- **Мигов Д.А.** Расчет надежности сети с ограничением на диаметр с применением точек сочленения // Автоматика и телемеханика. 2011. №7. С. 69–74.

See [docs/references/README.md](references/README.md) for full references and [docs/references/algorithm_checklist.md](references/algorithm_checklist.md) for implementation correspondence.

## Future Improvements

- Memoization/caching for repeated subproblems
- More efficient block decomposition algorithms
- GPU acceleration for parallel calculations
- Approximation algorithms for very large graphs
