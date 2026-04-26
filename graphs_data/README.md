# Graphs used in experiments

Graphs are stored in **KAO** format (compact CSR-style representation used
throughout the thesis literature). A full spec of the format lives in
[../docs/ALGORITHMS.md](../docs/ALGORITHMS.md); a one-line summary: a KAO
file has three comma-separated rows — `KAO` (cumulative degree prefix),
`FO` (flat 1-based adjacency list), and `Targets` (N marks indicating
which vertices are designated `s` and `t`).

Edge-list variants (`*_edgelist.txt`) exist for a few graphs as the
"source of truth" format used by the author of the original data set;
see [`Exceptions.h`](../include/graph_reliability/Exceptions.h) and
[`GraphOperations`](../src/GraphOperations.cpp) for the converter.

## Catalogue

### Complete graphs (smoke tests)

| File | Vertices | Edges | Purpose |
|---|---:|---:|---|
| `K4_kao.txt` | 4 | 6 | Minimal non-trivial test; every pair reachable in 1 hop |
| `K5_kao.txt` | 5 | 10 | Slightly larger smoke test |

### Synthetic "sausage" benchmarks

A sausage-k-MxM graph is k copies of the MxM mesh glued at a single
vertex each, producing a chain of blocks. The `s=0, t=N-1` pair lies at
opposite ends, so `dist(s, t) ≈ block_count × block_span` and the graph
stresses block-decomposition methods.

| File | Vertices | Edges | Blocks × Size |
|---|---:|---:|---|
| `2_3x3_blocks_kao.txt` | 17 | 24 | 2 × 3×3 |
| `3_blocks_sausage_3x3_kao.txt` | 29 | 48 | 3 × 3×3 |
| `3_blocks_sausage_4x4_kao.txt` | 50 | 88 | 3 × 4×4 (harder) |
| `4_blocks_sausage_3x3_kao.txt` | 39 | 66 | 4 × 3×3 |
| `5_blocks_sausage_3x3_kao.txt` | 49 | 84 | 5 × 3×3 |
| `6_blocks_sausage_3x3_kao.txt` | 59 | 102 | 6 × 3×3 |

### Real-world backbones

| File | Vertices | Edges | Source |
|---|---:|---:|---|
| `Geant2004_kao.txt` | 103 | 127 | GÉANT research-network backbone, 2004 snapshot |
| `Geant2009_kao.txt` | 390 | 503 | GÉANT backbone, 2009 snapshot (much larger) |
| `Geant2009_edgelist.txt` | 390 | 503 | Same as above in edge-list form |
| `IEEE-118-node_kao.txt` | 118 | 168 | IEEE 118-bus test system (power flow benchmark) |
| `UPS_of_Russia_composed_with_colored_cut_vertices_kao.txt` | 63 | 108 | Unified Power System of Russia (composed) |

### Working files

| File | Purpose |
|---|---|
| `input.txt` | Scratch pad for ad-hoc runs via `--run` |
| `GraphsToMerge.txt` / `MergedGraphs.txt` | Fixtures for the graph-composition utility |

## Choosing a test pair

Every KAO file has `Targets` marks: the first two `1`-positions are
treated as `s` and `t` by the auto-detection logic in
[`TestSuite::findSourceAndTargetVertices`](../src/TestSuite.cpp).
Override with explicit indices via `--run <file> <s> <t> <d> <method>`.

For `--cross-check`, diameter is **not** accepted manually — it is
auto-ranged from `dist(s, t)`; see
[../docs/EXPERIMENTS.md](../docs/EXPERIMENTS.md).
