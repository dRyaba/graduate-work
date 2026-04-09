/**
 * @file ReliabilityGraph.h
 * @brief Extended graph structure for reliability calculations with target vertices
 * @author Graduate Work Project
 * @date 2026
 * 
 * This class extends the base Graph with target vertices and specialized reliability
 * calculation methods. It provides four different algorithms for calculating network
 * reliability with diameter constraints:
 * 
 * - Method 0: Standard Factoring (baseline, no decomposition)
 * - Method 1: Recursive Decomposition (nested recursion, inefficient)
 * - Method 2: Simple Factoring with Decomposition (good balance)
 * - Method 3: M-Decomposition (fastest, recommended for production)
 * 
 * @see Graph for base graph functionality
 * @see docs/ALGORITHMS.md for detailed algorithm descriptions
 */

#pragma once

#include "Graph.h"
#include <optional>
#include <vector>
#include <utility>
#include <fstream>

namespace graph_reliability {

/**
 * @brief Result structure for reliability calculation methods
 */
struct ReliabilityResult {
    double reliability;           ///< Calculated reliability value
    long long recursions;         ///< Number of recursive calls
    double execution_time_sec;    ///< Execution time in seconds

    ReliabilityResult(double rel = 0.0, long long rec = 0, double time = 0.0)
        : reliability(rel), recursions(rec), execution_time_sec(time) {}
};

/**
 * @brief Extended graph structure for reliability calculations
 * 
 * This class extends the base Graph with target vertices and specialized
 * reliability calculation methods for network reliability analysis.
 */
class ReliabilityGraph : public Graph {
public:
    std::vector<int> target_vertices_; ///< Target vertices for reliability calculations

    /**
     * @brief Default constructor
     */
    ReliabilityGraph() = default;

    /**
     * @brief Constructor with CSR data and target vertices
     * @param kao Offset array for CSR format
     * @param fo Adjacency list
     * @param p_array Edge probabilities
     * @param targets Target vertices
     */
    ReliabilityGraph(std::vector<VertexId> kao,
                     std::vector<VertexId> fo,
                     std::vector<Probability> p_array,
                     std::vector<int> targets);

    /**
     * @brief Destructor
     */
    ~ReliabilityGraph() override = default;

    /**
     * @brief Calculate reliability with diameter constraint
     * @param diameter Maximum allowed diameter
     * @return Reliability result
     */
    ReliabilityResult calculateReliabilityWithDiameter(int diameter) const;

    /**
     * @brief Calculate reliability between two vertices with diameter constraint
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param diameter Maximum allowed diameter
     * @return Reliability result
     */
    ReliabilityResult calculateReliabilityBetweenVertices(VertexId source_vertex,
                                                         VertexId target_vertex,
                                                         int diameter) const;

    /**
     * @brief Calculate reliabilities for multiple diameters in ONE pass (modified factoring)
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param lower_bound Lower bound for diameter
     * @param upper_bound Upper bound for diameter
     * @return Pair of (cumulative reliabilities vector, recursion count)
     */
    std::pair<std::vector<double>, long long> calculateReliabilityMultipleDiameters(
        VertexId source_vertex,
        VertexId target_vertex,
        int lower_bound,
        int upper_bound) const;

    /**
     * @brief Calculate reliability with decomposition and simple factoring
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param upper_bound_diameter Upper bound for diameter
     * @return Reliability result
     */
    ReliabilityResult calculateReliabilityWithDecomposition(VertexId source_vertex,
                                                           VertexId target_vertex,
                                                           int upper_bound_diameter) const;

    /**
     * @brief Calculate reliability with M-decomposition (Method 3 - FASTEST)
     * @param source_vertex Source vertex (0-based)
     * @param target_vertex Target vertex (0-based)
     * @param upper_bound_diameter Upper bound for diameter constraint
     * @return Reliability result containing reliability value, recursion count, and execution time
     * 
     * This is the fastest method, using block decomposition with modified factoring.
     * It computes reliabilities for multiple diameters simultaneously, avoiding
     * redundant calculations.
     * 
     * @complexity O(B × E × D) where B=blocks, E=edges, D=diameter range
     * @note Recommended for production use
     * @note Uses block decomposition to reduce problem size
     * @see calculateReliabilityMultipleDiameters for the core optimization
     */
    ReliabilityResult calculateReliabilityWithMDecomposition(VertexId source_vertex,
                                                            VertexId target_vertex,
                                                            int upper_bound_diameter) const;

    /**
     * @brief Calculate reliability with parallel M-decomposition
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param upper_bound_diameter Upper bound for diameter
     * @return Reliability result
     */
    ReliabilityResult calculateReliabilityWithParallelMDecomposition(VertexId source_vertex,
                                                                    VertexId target_vertex,
                                                                    int upper_bound_diameter) const;

    /**
     * @brief Calculate reliability with recursive decomposition
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param upper_bound_diameter Upper bound for diameter
     * @return Reliability result
     */
    ReliabilityResult calculateReliabilityWithRecursiveDecomposition(VertexId source_vertex,
                                                                    VertexId target_vertex,
                                                                    int upper_bound_diameter) const;

    /**
     * @brief Calculate reliability with Cancela-Petingi method (Method 4 - path-based factoring)
     * @param source_vertex Source vertex (0-based)
     * @param target_vertex Target vertex (0-based)
     * @param diameter Maximum path length (diameter constraint)
     * @return Reliability result
     *
     * Uses path lists instead of graphs; supports ESS and ISPT for acceleration.
     * @see TR0103.pdf (Cancela, Petingi), NesterovMigov.pdf
     */
    ReliabilityResult calculateReliabilityCancelaPetingi(VertexId source_vertex,
                                                        VertexId target_vertex,
                                                        int diameter) const;

    /**
     * @brief Calculate reliability CDF with Cancela-Petingi for multiple diameters
     * @param source_vertex Source vertex (0-based)
     * @param target_vertex Target vertex (0-based)
     * @param lower_bound Lower bound for diameter
     * @param upper_bound Upper bound for diameter
     * @return Pair of (cumulative reliabilities vector, recursion count)
     *
     * Multi-diameter CPFM: computes R(d) for d in [lower_bound, upper_bound] in ONE pass.
     * Used in M-Decomposition + CPFM hybrid (Method 5).
     */
    std::pair<std::vector<double>, long long> calculateReliabilityCancelaPetingiMulti(
        VertexId source_vertex,
        VertexId target_vertex,
        int lower_bound,
        int upper_bound) const;

    /**
     * @brief Calculate reliability with M-Decomposition + CPFM hybrid (Method 5)
     * @param source_vertex Source vertex (0-based)
     * @param target_vertex Target vertex (0-based)
     * @param upper_bound_diameter Upper bound for diameter constraint
     * @return Reliability result
     *
     * Combines block decomposition (M-Decomposition) with path-based factoring (CPFM).
     * Uses CPFM inside each block for efficiency, then convolves block reliabilities.
     *
     * @complexity O(B × 2^E_block × paths) where B=blocks, E_block=edges per block
     * @see calculateReliabilityCancelaPetingiMulti for the core CPFM optimization
     */
    ReliabilityResult calculateReliabilityWithMDecompositionCPFM(VertexId source_vertex,
                                                                 VertexId target_vertex,
                                                                 int upper_bound_diameter) const;

    /**
     * @brief Check if graph is k-connected
     * @return True if graph is k-connected, false otherwise
     */
    bool isKConnected() const;

    /**
     * @brief Check distance constraint using Floyd-Warshall algorithm
     * @param max_distance Maximum allowed distance
     * @return True if distance constraint is satisfied
     */
    bool checkDistanceConstraint(int max_distance) const;

    /**
     * @brief Decompose graph into blocks
     * @return Block decomposition result
     */
    std::vector<int> decomposeIntoBlocks() const;

    /**
     * @brief Decompose graph into k-blocks
     * @return K-block decomposition result
     */
    std::vector<int> decomposeIntoKBlocks() const;

    /**
     * @brief Restore block from decomposition
     * @param block_number Block number to restore
     * @param decomposition Decomposition result
     * @return Restored block graph
     */
    ReliabilityGraph restoreBlock(int block_number, const std::vector<int>& decomposition) const;

    /**
     * @brief Restore k-block from decomposition
     * @param block_number Block number to restore
     * @param decomposition Decomposition result
     * @return Restored k-block graph
     */
    ReliabilityGraph restoreKBlock(int block_number, const std::vector<int>& decomposition) const;

    /**
     * @brief Check connectivity without specified block
     * @param block_number Block number to exclude
     * @param block_decomposition Block decomposition vector
     * @return True if connected without the block
     */
    bool checkConnectivityWithoutBlock(int block_number, const std::vector<int>& block_decomposition) const;

    /**
     * @brief Remove edge from reliability graph
     * @param from Source vertex
     * @param to Destination vertex
     * @return New reliability graph without the edge
     */
    ReliabilityGraph removeEdge(VertexId from, VertexId to) const;

    /**
     * @brief Create induced subgraph
     * @param vertices Vertices to include in subgraph
     * @param block_number Block number for filtering
     * @return Induced subgraph
     */
    ReliabilityGraph createInducedSubgraph(const std::vector<VertexId>& vertices, int block_number) const;

    /**
     * @brief Find last unreliable edge in the graph
     * @return Index of last unreliable edge or nullopt if all edges are reliable
     */
    std::optional<EdgeId> findLastUnreliableEdge() const;

    /**
     * @brief Find reverse edge indices for a given edge
     * @param edge_index Index of the edge
     * @return Pair of (source_vertex, reverse_edge_index)
     */
    std::pair<VertexId, EdgeId> findReverseEdgeIndices(EdgeId edge_index) const;

    /**
     * @brief Check if vertices are in the same block
     * @param vertex_set Set of vertices
     * @param block_number Block number
     * @param vertex1 First vertex
     * @param vertex2 Second vertex
     * @return True if vertices are in the same block
     */
    bool areVerticesInSameBlock(const std::vector<VertexId>& vertex_set,
                               int block_number,
                               VertexId vertex1,
                               VertexId vertex2) const;

    /**
     * @brief Output graph to file stream
     * @param output_stream Output file stream
     */
    void outputToFile(std::ofstream& output_stream) const;

    /**
     * @brief Swap vertices in the graph
     * @param vertex1 First vertex
     * @param vertex2 Second vertex
     * @return New graph with swapped vertices
     */
    ReliabilityGraph swapVertices(VertexId vertex1, VertexId vertex2) const;

    /**
     * @brief Extract block graph and vertex mappings
     * @param block_id Block ID
     * @param decomposition Block decomposition vector
     * @param out_graph Output block subgraph
     * @param map_new_to_orig Mapping from new to original vertex IDs
     * @param map_orig_to_new Mapping from original to new vertex IDs
     */
    void getBlockGraphAndMap(
        int block_id,
        const std::vector<int>& decomposition,
        ReliabilityGraph& out_graph,
        std::vector<int>& map_new_to_orig,
        std::vector<int>& map_orig_to_new
    ) const;

    /**
     * @brief Get list of blocks containing a specific vertex
     * @param v Vertex ID
     * @param decomposition Block decomposition vector
     * @return List of block IDs
     */
    std::vector<int> getBlocksContainingVertex(VertexId v, const std::vector<int>& decomposition) const;

    /**
     * @brief Calculate minimum path length (d_k^0) within a block
     * @param block_id Block ID
     * @param decomposition Block decomposition vector
     * @param entry Entry vertex (original graph index)
     * @param exit Exit vertex (original graph index)
     * @return Minimum path length from entry to exit within the block, or -1 if unreachable
     * 
     * Uses BFS to find shortest path in the block graph.
     * Used in gap optimization for Method 5.
     */
    int calculateMinBlockDiameter(int block_id, const std::vector<int>& decomposition,
                                   VertexId entry, VertexId exit) const;

private:
    /**
     * @brief Result of block chain decomposition for a given (s, t) pair.
     *
     * Contains everything downstream methods need:
     *   - decomposition   : block labels from decomposeIntoBlocks()
     *   - ordered_block_ids: block IDs along the s→t path
     *   - final_aps       : entry/exit vertices per block,
     *                       final_aps[0]=s, final_aps[i]=AP between block i-1 and i,
     *                       final_aps[N]=t
     *   - num_blocks      : total number of blocks in the graph
     *   - valid           : false when single-block graph, s/t unreachable, or path not found
     */
    struct BlockChain {
        std::vector<int> decomposition;
        std::vector<int> ordered_block_ids;
        std::vector<int> final_aps;
        int num_blocks = 0;
        bool valid     = false;
    };

    /** @brief Build block chain for a given (s, t) pair. */
    BlockChain buildBlockChain(VertexId s, VertexId t) const;

    /**
     * @brief Recursive solver for block chain reliability (uses MODIFIED factoring)
     */
    std::vector<double> solveRecursiveForBlockChain(
        int current_idx,
        int entry_node_orig,
        int target_node_orig,
        int max_len,
        const std::vector<int>& decomposition,
        const std::vector<int>& block_ids,
        const std::vector<int>& aps_orig,
        long long& recursion_counter
    ) const;

    /**
     * @brief TRUE NESTED RECURSION solver for Level 1 (Recursive Decomposition)
     * 
     * This implements the theoretically correct but INEFFICIENT nested recursion.
     * For each path length L in block_i, makes a RECURSIVE CALL to block_{i+1}.
     * Returns a SCALAR reliability value (not a vector).
     */
    double solveNestedRecursive(
        int current_idx,
        int entry_node_orig,
        int target_node_orig,
        int remaining_diameter,
        const std::vector<int>& decomposition,
        const std::vector<int>& block_ids,
        const std::vector<int>& aps_orig,
        long long& recursion_counter
    ) const;

    /**
     * @brief Internal factoring algorithm
     */
    void performFactoring(ReliabilityGraph graph, int variant, int diameter, double reliability) const;

    /**
     * @brief Internal 2-vertex factoring algorithm
     */
    void perform2VertexFactoring(ReliabilityGraph graph,
                                VertexId source_vertex,
                                VertexId target_vertex,
                                int variant,
                                int diameter,
                                double reliability) const;

    /**
     * @brief Internal M-factoring algorithm
     */
    void performMFactoring(ReliabilityGraph graph,
                          VertexId source_vertex,
                          VertexId target_vertex,
                          int variant,
                          int diameter,
                          double reliability,
                          int lower_bound,
                          int upper_bound) const;

    /**
     * @brief Internal parallel M-factoring algorithm
     */
    void performParallelMFactoring(ReliabilityGraph graph,
                                  VertexId source_vertex,
                                  VertexId target_vertex,
                                  int variant,
                                  int diameter,
                                  double reliability,
                                  int lower_bound,
                                  int upper_bound) const;
};

} // namespace graph_reliability

