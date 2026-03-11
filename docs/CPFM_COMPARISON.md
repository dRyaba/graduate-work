# CPFM Implementation: Cancela-Petingi Factoring Method

This document describes the optimized Cancela-Petingi implementation in graduate-work.

## Bug Fix: ISPT Comparison

The original ISPT (Improved Series-Parallel Transform) had a bug: it compared P(e) sets using **all paths**, but should only compare **feasible paths**.

### Problem

When comparing P(e) and P(f) to decide if edges can be merged:

- Original: compared all path indices in P(e) and P(f)
- Correct: compare only feasible paths (paths not yet eliminated by delete operations)

This caused incorrect reliability values on larger graphs.

### Solution

Modified `applyISPT()` and `applyGlobalISPT()` to filter by `feasible_[pi]` before comparing.

## Final Results

All values now match the reference implementation (scientific advisor's program).

### K4 (s=0, t=3, D=3, p=0.9)

| Metric | Value |
|--------|-------|
| Reliability | 0.997848000000000 |
| Recursions | 17 |

### K5 (s=0, t=4, D=4, p=0.9)

| Metric | Value |
|--------|-------|
| Reliability | 0.999794802600000 |
| Recursions | 206 |

### 2_3x3_blocks (s=0, t=16, D=10, p=0.9)

| Metric | Value |
|--------|-------|
| Reliability | 0.945586832221472 |
| Recursions | 5688 |

### 3_blocks_sausage_3x3 (s=0, t=28, D=12, p=0.9)

| Metric | Value |
|--------|-------|
| Reliability | **0.934403607310576** |
| Recursions | 15221 |

**Reference value (scientific advisor):** 0.934403607310576 — **exact match**

## Performance Comparison

| Graph | Without ISPT | With ISPT | Improvement |
|-------|--------------|-----------|-------------|
| K4 | 21 | 17 | -19% |
| K5 | 254 | 206 | -19% |
| 2_3x3_blocks | 31847 | 5688 | **-82%** |
| 3_blocks_sausage_3x3 | 66981 | 15221 | **-77%** |

## Algorithm Components

### 1. Path Enumeration

`PathEnumerator::enumeratePaths()` — DFS enumeration of all simple paths from s to t with length ≤ D.

### 2. Global ISPT (before recursion)

`applyGlobalISPT()` — Merge all edges with identical feasible P(e) sets. This is more aggressive than per-pivot ISPT.

### 3. Dynamic Pivot Selection (ESS)

`selectPivotEdge()` — Select edge with maximum |P(e)| among feasible paths.

### 4. Per-Pivot ISPT

`applyISPT(e)` — After selecting pivot e, merge additional edges with same feasible P(e).

### 5. Factoring

```
facto(state):
    if connected: return 1.0
    if no feasible paths: return 0.0
    e = selectPivot()
    applyISPT(e)
    re = reliability(e)
    return re * facto(contract(e)) + (1-re) * facto(delete(e))
```

## Test Commands

```powershell
# Build
cmake --build build --target GraphReliabilityAnalysis

# Test with reference values
.\build\graph_reliability.exe --run graphs_data\3_blocks_sausage_3x3_kao.txt 0 28 12 4 1
# Expected: 0.934403607310576

.\build\graph_reliability.exe --run graphs_data\K4_kao.txt 0 3 3 4 1
# Expected: 0.997848

.\build\graph_reliability.exe --run graphs_data\K5_kao.txt 0 4 4 4 1
# Expected: 0.9997948026
```

## Files Modified

- `src/CancelaPetingiState.cpp` — Fixed ISPT to compare feasible paths only
- `src/ReliabilityGraph.cpp` — CPFM with global ISPT + per-pivot ISPT
- `tests/test_ReliabilityGraph.cpp` — Updated expected recursion counts
