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

### Method 4: Cancela-Petingi (Level 4)

**Description**: Path-based factoring method. Operates on path lists instead of graphs; uses series-parallel transformation (SPT), edge selection strategy (ESS), and improved SPT on each recursive call (ISPT).

**Algorithm**:
1. Apply series transformation: remove 2-degree non-terminal chains
2. Enumerate all paths Pst(d) of length ≤ d between s and t
3. Build P(e) for each edge: paths containing e
4. FACTO: select pivot e (max |P(e)|), apply ISPT, Contract/Delete branches
5. R = re × RContract + (1-re) × RDelete

**Key Optimizations**:
- ESS: pivot edge with maximum |P(e)|
- SPT: chain reduction before path enumeration
- ISPT: merge edges with P(e)=P(f) on each recursive call

**Complexity**: O(2^m × paths) with early termination; typically fewer recursions than graph-based factoring

**Use Case**: Alternative to graph-based methods; single-block calculations

### Method 5: M-Decomposition + CPFM (Level 5)

**Description**: Hybrid method combining block decomposition (M-Decomposition) with path-based factoring (Cancela-Petingi). Uses CPFM inside each block for efficiency, then convolves block reliabilities using Migov's formula.

**Algorithm**:
1. Decompose graph into k-blocks
2. Build block graph and find path through it
3. For each block in path:
   - Enumerate all paths within the block
   - Apply global ISPT for series-parallel transformation
   - Use multi-diameter CPFM to compute CDF for all diameters in one pass
   - CDF[d] = P(exists path of length ≤ d)
4. Convert CDF to PMF and convolve block reliabilities

**Key Innovations**:
- Multi-diameter CPFM: computes reliability for all diameters [0, D] simultaneously
- No early termination on success: ensures correct CDF computation for all lengths
- ISPT acceleration: reduces number of edges to factor

**Complexity**: O(B × 2^E_block × paths) where B=blocks, E_block=edges per block

**Advantages**:
- Fewer recursions than Method 3 due to ISPT
- Exact computation (no approximation)
- Combines benefits of decomposition and path-based approach

**Use Case**: Block-structured graphs, production use with ISPT benefits

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
| Cancela-Petingi | O(2^m × paths) | Single-block graphs |
| M-Decomp + CPFM | O(B × 2^E × paths) | Block-structured with ISPT |

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
- **Cancela H., Petingi L.** Diameter constrained network reliability: exact evaluation by factorization and bounds. (Path-based formulation, recursive factoring with path lists)
- **Nesterov S., Migov D.** Series-parallel transformation for diameter constrained network reliability computation. (SPT, ESS, ISPT on each recursive call of CPFM)

See [docs/references/README.md](references/README.md) for full references and [docs/references/algorithm_checklist.md](references/algorithm_checklist.md) for implementation correspondence.

## Experimental Results

All benchmarks were run on a single CPU core (Release build, -O3). Timeout = 300 s.  
"d_min" denotes the minimum diameter for which a path exists between the chosen s–t pair.

### Test Graph Properties

| Graph | V | E | Blocks | Bridges | Max block edges |
|---|---|---|---|---|---|
| 2-block 3×3 sausage | 18 | 26 | 2 | 0 | 12 |
| 3-block 3×3 sausage | 27 | 39 | 3 | 0 | 12 |
| 3-block 4×4 sausage | 38 | 60 | 3 | 0 | 24 |
| 4-block 3×3 sausage | 36 | 52 | 4 | 0 | 12 |
| 5-block 3×3 sausage | 45 | 65 | 5 | 0 | 12 |
| 6-block 3×3 sausage | 54 | 78 | 6 | 0 | 12 |
| GEANT 2004 | 103 | 127 | 50 | 41 | 25 |
| GEANT 2009 | 390 | 503 | 91 | 75 | 238 |
| IEEE 118-node | 118 | 168 | 14 | 12 | 144 |
| UPS Russia | 63 | 108 | 17 | 12 | 85 |

### Sausage Chain Graphs — at d = d_min (gap = 0)

| Graph | d | m3 | m4 | **m5** | Speedup m5/m3 |
|---|---|---|---|---|---|
| 2-block 3×3 | 8 | 16 ms | 37 ms | **2.6 ms** | ×6 |
| 3-block 3×3 | 10 | 565 ms | 0.08 ms | **0.11 ms** | **×4961** |
| 3-block 4×4 | 13 | TIMEOUT | 0.86 ms | **0.86 ms** | — |
| 4-block 3×3 | 14 | 1137 ms | 0.82 ms | **0.17 ms** | **×6741** |
| 5-block 3×3 | 18 | 1562 ms | 12.7 ms | **0.20 ms** | **×7854** |
| 6-block 3×3 | 22 | 2056 ms | 208 ms | **0.26 ms** | **×7880** |

### Sausage Chain Graphs — effect of increasing d (gap > 0)

Key finding: m4 degrades exponentially for d > d_min on sausage chains. m5 degrades gracefully.

| Graph | d | gap | m3 | m4 | m5 | m5/m3 |
|---|---|---|---|---|---|---|
| 3-block 3×3 | 10 | 0 | 565 ms | 0.08 ms | **0.11 ms** | ×4961 |
| 3-block 3×3 | 11 | 1 | 564 ms | 11 ms | 305 ms | ×2 |
| 3-block 3×3 | 12 | 2 | 562 ms | 763 ms | 450 ms | ×1.3 |
| 5-block 3×3 | 18 | 0 | 1562 ms | 12.7 ms | **0.20 ms** | ×7854 |
| 5-block 3×3 | 19 | 1 | 1563 ms | **TIMEOUT** | 904 ms | ×1.7 |
| 5-block 3×3 | 20 | 2 | 1569 ms | **TIMEOUT** | 1278 ms | ×1.2 |
| 6-block 3×3 | 22 | 0 | 2056 ms | 208 ms | **0.26 ms** | ×7880 |
| 6-block 3×3 | 23 | 1 | 2067 ms | **TIMEOUT** | 1189 ms | ×1.7 |
| 6-block 3×3 | 24 | 2 | 2058 ms | **TIMEOUT** | 1688 ms | ×1.2 |

All methods produce identical reliability values (tolerance < 1e-9).

**Observations:**
- At gap=0 (d = d_min): m5 is ×4961–×7880 faster than m3. Gap optimisation reduces per-block CPFM to a single CDF point.
- At gap=1: m4 begins to time out on longer chains; m5 is ~×2 vs m3.
- At gap=2: m4 completely unavailable for chains ≥ 5 blocks; m5 ≈ m3 (gap budget spread across all blocks makes per-block CPFM expensive).
- **m5 is the only method applicable for all d values on sausage chains.** It is never worse than m3, and vastly superior at tight diameter constraints.

### Real-World Networks

| Graph | d | m3 | **m4** | m5 | Note |
|---|---|---|---|---|---|
| GEANT 2004 | 6 | 3 ms | **0.03 ms** | 0.15 ms | Small blocks (max 25 e), no fallback |
| GEANT 2004 | 8 | — | 0.37 ms | 1.55 ms | Per-block CPFM costlier than global |
| GEANT 2009 | 12 | TIMEOUT | **0.04 ms** | 0.04 ms | Fallback to m4 (block 238 e) |
| IEEE 118-node | 10 | 36.8 s | **0.04 ms** | 0.04 ms | Fallback to m4 (block 144 e) |
| IEEE 118-node | 12 | — | 8.8 ms | 9.3 ms | — |
| UPS Russia | 18 | TIMEOUT | 158 ms | 160 ms | Fallback to m4 (block 85 e) |

**Observations:**
- For graphs with large biconnected components (IEEE-118, UPS, GEANT2009): m5 auto-falls back to m4; both methods behave identically.
- GEANT2004 has small blocks (max 25 edges): m4 wins here because global ESS/ISPT is more effective than per-block CPFM at d > d_min.
- m3 is completely impractical on real networks with large blocks.

### Conclusion

| Use case | Recommended method | Reason |
|---|---|---|
| Sausage chains, d = d_min | **m5** | Gap=0 optimisation: ×4961–×7880 vs m3 |
| Sausage chains, d > d_min | **m5** | m4 times out; m5 degrades gracefully |
| Real networks, large blocks | **m4** | Global ESS/ISPT more effective per-block CPFM |
| Real networks, small blocks | **m4** | Global CPFM beats per-block at d > d_min |
| Universal safe choice | **m5** | Automatically falls back to m4 when blocks are large |

| Use case | Recommended method |
|---|---|
| Sausage chains / multi-block structured graphs | **m5** |
| Real networks with large biconnected components | **m4** |
| Small single-block graphs | m4 (or m0 for reference) |
| Theoretical comparison | m3 (modified factoring baseline) |

## Future Improvements

- Memoization/caching for repeated subproblems
- More efficient block decomposition algorithms
- GPU acceleration for parallel calculations
- Approximation algorithms for very large graphs
