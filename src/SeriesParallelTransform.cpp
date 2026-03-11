/**
 * @file SeriesParallelTransform.cpp
 * @brief Series transformation for DCR (chain removal)
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/SeriesParallelTransform.h"
#include "graph_reliability/Graph.h"
#include <map>
#include <set>
#include <vector>
#include <algorithm>

namespace graph_reliability {

namespace {

struct EdgeKey {
    int u, v;
    EdgeKey(int a, int b) : u(std::min(a, b)), v(std::max(a, b)) {}
    bool operator<(const EdgeKey& o) const {
        return u < o.u || (u == o.u && v < o.v);
    }
};

using EdgeMap = std::map<EdgeKey, double>;

Graph buildGraphFromEdges(const EdgeMap& edges, int num_vertices) {
    std::vector<std::vector<std::pair<int, double>>> adj(num_vertices);
    for (const auto& [key, r] : edges) {
        adj[key.u].emplace_back(key.v, r);
        adj[key.v].emplace_back(key.u, r);
    }

    std::vector<Graph::VertexId> kao, fo;
    std::vector<Graph::Probability> p_array;
    kao.push_back(0);

    for (int u = 0; u < num_vertices; ++u) {
        for (const auto& [v, r] : adj[u]) {
            fo.push_back(v);
            p_array.push_back(r);
        }
        kao.push_back(static_cast<Graph::VertexId>(fo.size()));
    }

    return Graph(kao, fo, p_array);
}

EdgeMap graphToEdgeMap(const Graph& g) {
    EdgeMap em;
    size_t n = g.numVertices();
    for (size_t u = 0; u < n; ++u) {
        for (size_t k = g.kao_[u]; k < g.kao_[u + 1]; ++k) {
            int v = g.fo_[k];
            if (static_cast<int>(u) < v) {
                em[EdgeKey(u, v)] = g.p_array_[k];
            }
        }
    }
    return em;
}

} // namespace

std::pair<Graph, double> applySeriesTransformation(const Graph& graph,
                                                   const std::vector<int>& terminals) {
    size_t n = graph.numVertices();
    std::set<int> term_set(terminals.begin(), terminals.end());

    EdgeMap edges = graphToEdgeMap(graph);
    double factor = 1.0;
    bool changed = true;

    while (changed) {
        changed = false;
        for (int v = 0; v < static_cast<int>(n); ++v) {
            if (term_set.count(v)) continue;

            std::vector<std::pair<int, double>> neighbors;
            for (const auto& [key, r] : edges) {
                if (key.u == v) neighbors.emplace_back(key.v, r);
                else if (key.v == v) neighbors.emplace_back(key.u, r);
            }
            if (neighbors.size() != 2) continue;

            int s = neighbors[0].first;
            int t = neighbors[1].first;
            double r_sv = neighbors[0].second;
            double r_vt = neighbors[1].second;

            EdgeKey sv_key(s, v);
            EdgeKey vt_key(v, t);
            EdgeKey st_key(s, t);

            if (edges.count(st_key)) continue;

            edges.erase(sv_key);
            edges.erase(vt_key);
            edges[st_key] = r_sv * r_vt;
            changed = true;
            break;
        }
    }

    return {buildGraphFromEdges(edges, static_cast<int>(n)), factor};
}

} // namespace graph_reliability
