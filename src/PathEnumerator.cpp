/**
 * @file PathEnumerator.cpp
 * @brief Implementation of path enumeration for Cancela-Petingi method
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/PathEnumerator.h"
#include <algorithm>

namespace graph_reliability {

PathEnumerator::EdgeId PathEnumerator::getCanonicalEdgeId(const Graph& graph, VertexId u, VertexId v) {
    VertexId a = std::min(u, v);
    VertexId b = std::max(u, v);
    return static_cast<EdgeId>(graph.findEdge(a, b));
}

void PathEnumerator::dfsEnumerate(const Graph& graph,
                                  VertexId current,
                                  VertexId target,
                                  int max_diameter,
                                  int current_length,
                                  std::vector<bool>& visited,
                                  std::vector<EdgeId>& current_path_edges,
                                  std::vector<Path>& result) {
    if (current_length > max_diameter) return;

    for (size_t i = graph.kao_[current]; i < graph.kao_[current + 1]; ++i) {
        VertexId next = graph.fo_[i];
        EdgeId edge_id = getCanonicalEdgeId(graph, current, next);

        if (next == target) {
            if (current_length + 1 <= max_diameter) {
                current_path_edges.push_back(edge_id);
                result.push_back(current_path_edges);
                current_path_edges.pop_back();
            }
            continue;
        }

        if (!visited[next] && current_length + 1 < max_diameter) {
            visited[next] = true;
            current_path_edges.push_back(edge_id);
            dfsEnumerate(graph, next, target, max_diameter, current_length + 1,
                         visited, current_path_edges, result);
            current_path_edges.pop_back();
            visited[next] = false;
        }
    }
}

std::vector<PathEnumerator::Path> PathEnumerator::enumeratePaths(const Graph& graph,
                                                                VertexId source,
                                                                VertexId target,
                                                                int max_diameter) {
    std::vector<Path> result;

    if (max_diameter < 1) return result;

    size_t n = graph.numVertices();
    if (source < 0 || source >= static_cast<VertexId>(n) ||
        target < 0 || target >= static_cast<VertexId>(n)) {
        return result;
    }

    if (source == target) {
        return result;
    }

    std::vector<bool> visited(n, false);
    visited[source] = true;
    std::vector<EdgeId> current_path_edges;

    dfsEnumerate(graph, source, target, max_diameter, 0,
                 visited, current_path_edges, result);

    return result;
}

} // namespace graph_reliability
