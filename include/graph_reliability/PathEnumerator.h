/**
 * @file PathEnumerator.h
 * @brief Path enumeration for Cancela-Petingi method
 * @author Graduate Work Project
 * @date 2026
 *
 * Enumerates all simple paths of length at most d between source and target.
 * Used by CPFM (Cancela-Petingi Factoring Method) for path-based reliability.
 *
 * @see ReliabilityGraph::calculateReliabilityCancelaPetingi
 * @see TR0103.pdf (Cancela, Petingi)
 */

#pragma once

#include "Graph.h"
#include <vector>

namespace graph_reliability {

/**
 * @brief Enumerates paths between two vertices for diameter-constrained reliability
 *
 * Uses DFS to find all simple paths of length at most max_diameter.
 * Each path is represented as a sequence of canonical edge indices.
 */
class PathEnumerator {
public:
    using EdgeId = Graph::EdgeId;
    using VertexId = Graph::VertexId;

    /**
     * @brief Path as sequence of canonical edge indices
     * Canonical edge (u,v) with u < v uses findEdge(u,v)
     */
    using Path = std::vector<EdgeId>;

    /**
     * @brief Enumerate all simple paths of length <= max_diameter between source and target
     * @param graph Input graph (CSR format)
     * @param source Source vertex (0-based)
     * @param target Target vertex (0-based)
     * @param max_diameter Maximum path length (number of edges)
     * @return List of paths; each path is a vector of canonical edge indices
     *
     * @complexity O(delta^d) where delta is max degree, d is max_diameter
     */
    static std::vector<Path> enumeratePaths(const Graph& graph,
                                            VertexId source,
                                            VertexId target,
                                            int max_diameter);

    /**
     * @brief Get canonical edge index for undirected edge (u, v)
     * Uses the index from the smaller vertex for consistency
     */
    static EdgeId getCanonicalEdgeId(const Graph& graph, VertexId u, VertexId v);

private:
    static void dfsEnumerate(const Graph& graph,
                             VertexId current,
                             VertexId target,
                             int max_diameter,
                             int current_length,
                             std::vector<bool>& visited,
                             std::vector<EdgeId>& current_path_edges,
                             std::vector<Path>& result);
};

} // namespace graph_reliability
