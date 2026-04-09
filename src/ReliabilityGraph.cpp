/**
 * @file ReliabilityGraph.cpp
 * @brief Implementation of the ReliabilityGraph class
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/ReliabilityGraph.h"
#include "graph_reliability/Logger.h"
#include "graph_reliability/PathEnumerator.h"
#include "graph_reliability/CancelaPetingiState.h"
#include "graph_reliability/SeriesParallelTransform.h"
#include <algorithm>
#include <functional>
#include <queue>
#include <stack>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <map>
#include <set>
#include <cmath>
#include <chrono>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace graph_reliability {

// Profiling structure for Method 5
struct Method5Profile {
    double decomposition_time = 0.0;
    double block_graph_build_time = 0.0;
    double path_finding_time = 0.0;
    double path_enumeration_time = 0.0;
    double global_ispt_time = 0.0;
    double factoring_time = 0.0;
    double convolution_time = 0.0;
    double total_time = 0.0;
    long long total_recursions = 0;
    long long total_paths = 0;
    int num_blocks = 0;
    std::vector<int> paths_per_block;
    std::vector<double> time_per_block;
    
    void print() const {
        std::cerr << "=== Method 5 Profile ===" << std::endl;
        std::cerr << "Total time: " << std::fixed << std::setprecision(6) << total_time << "s" << std::endl;
        if (total_time > 0) {
            std::cerr << "  Decomposition:      " << std::setw(10) << decomposition_time << "s (" 
                      << std::setprecision(1) << 100.0 * decomposition_time / total_time << "%)" << std::endl;
            std::cerr << "  Block graph build:  " << std::setw(10) << std::setprecision(6) << block_graph_build_time << "s (" 
                      << std::setprecision(1) << 100.0 * block_graph_build_time / total_time << "%)" << std::endl;
            std::cerr << "  Path finding:       " << std::setw(10) << std::setprecision(6) << path_finding_time << "s (" 
                      << std::setprecision(1) << 100.0 * path_finding_time / total_time << "%)" << std::endl;
            std::cerr << "  Path enumeration:   " << std::setw(10) << std::setprecision(6) << path_enumeration_time << "s (" 
                      << std::setprecision(1) << 100.0 * path_enumeration_time / total_time << "%)" << std::endl;
            std::cerr << "  Global ISPT:        " << std::setw(10) << std::setprecision(6) << global_ispt_time << "s (" 
                      << std::setprecision(1) << 100.0 * global_ispt_time / total_time << "%)" << std::endl;
            std::cerr << "  Factoring:          " << std::setw(10) << std::setprecision(6) << factoring_time << "s (" 
                      << std::setprecision(1) << 100.0 * factoring_time / total_time << "%)" << std::endl;
            std::cerr << "  Convolution:        " << std::setw(10) << std::setprecision(6) << convolution_time << "s (" 
                      << std::setprecision(1) << 100.0 * convolution_time / total_time << "%)" << std::endl;
        }
        std::cerr << "Num blocks: " << num_blocks << std::endl;
        std::cerr << "Total recursions: " << total_recursions << std::endl;
        std::cerr << "Total paths: " << total_paths << std::endl;
        for (size_t i = 0; i < paths_per_block.size(); ++i) {
            std::cerr << "  Block " << i << ": " << paths_per_block[i] << " paths, " 
                      << std::setprecision(6) << time_per_block[i] << "s" << std::endl;
        }
        std::cerr << "========================" << std::endl;
    }
};

// Global profiling flag - enable with -DPROFILE_METHOD5
#ifdef PROFILE_METHOD5
static thread_local Method5Profile g_profile;
static thread_local bool g_profiling_enabled = true;
#define PROFILE_START(var) auto var##_start = std::chrono::high_resolution_clock::now()
#define PROFILE_END(var, field) if (g_profiling_enabled) { \
    auto var##_end = std::chrono::high_resolution_clock::now(); \
    g_profile.field += std::chrono::duration<double>(var##_end - var##_start).count(); \
}
#define PROFILE_ADD_PATHS(n) if (g_profiling_enabled) g_profile.total_paths += (n)
#define PROFILE_PRINT() if (g_profiling_enabled) g_profile.print()
#define PROFILE_RESET() g_profile = Method5Profile{}
#else
#define PROFILE_START(var) ((void)0)
#define PROFILE_END(var, field) ((void)0)
#define PROFILE_ADD_PATHS(n) ((void)0)
#define PROFILE_PRINT() ((void)0)
#define PROFILE_RESET() ((void)0)
#endif

// -----------------------------------------------------------------------------
// Constructor & Basic Methods
// -----------------------------------------------------------------------------

ReliabilityGraph::ReliabilityGraph(std::vector<VertexId> kao,
                                   std::vector<VertexId> fo,
                                   std::vector<Probability> p_array,
                                   std::vector<int> targets)
    : Graph(std::move(kao), std::move(fo), std::move(p_array)),
      target_vertices_(std::move(targets)) {
}

std::optional<Graph::EdgeId> ReliabilityGraph::findLastUnreliableEdge() const {
    for (int i = static_cast<int>(p_array_.size()) - 1; i >= 0; --i) {
        if (p_array_[i] < 1.0) return i;
    }
    return std::nullopt;
}

std::pair<Graph::VertexId, Graph::EdgeId> ReliabilityGraph::findReverseEdgeIndices(EdgeId edge_index) const {
    VertexId from = 0;
    size_t n = numVertices();
    // Find 'from' vertex such that edge_index belongs to its adjacency list
    // Binary search is possible since KAO is sorted, but linear is fine for now
    while (from < n && kao_[from + 1] <= edge_index) {
        from++;
    }
    
    VertexId to = fo_[edge_index];
    // Search in 'to's edges for 'from'
    for (size_t k = kao_[to]; k < kao_[to + 1]; ++k) {
        if (fo_[k] == from) return {from, static_cast<EdgeId>(k)};
    }
    
    return {from, -1};
}

ReliabilityGraph ReliabilityGraph::removeEdge(VertexId from, VertexId to) const {
    Graph g = Graph::removeEdge(from, to);
    return ReliabilityGraph(g.kao_, g.fo_, g.p_array_, target_vertices_);
}

ReliabilityGraph ReliabilityGraph::swapVertices(VertexId vertex1, VertexId vertex2) const {
    Graph g = Graph::swapVertices(vertex1, vertex2);
    
    std::vector<int> new_targets = target_vertices_;
    size_t n = numVertices();
    if (vertex1 < n && vertex2 < n) {
        std::swap(new_targets[vertex1], new_targets[vertex2]);
    }
    
    return ReliabilityGraph(g.kao_, g.fo_, g.p_array_, new_targets);
}

bool ReliabilityGraph::isKConnected() const {
    int target_count = 0;
    VertexId start_node = -1;
    size_t n = numVertices();
    
    for (size_t i = 0; i < n; ++i) {
        if (target_vertices_[i] == 1) {
            target_count++;
            if (start_node == -1) start_node = i;
        }
    }

    if (target_count <= 1) return true;
    if (start_node == -1) return false;

    std::vector<int> spot(n, 0); 
    std::queue<VertexId> q;
    
    q.push(start_node);
    spot[start_node] = 1;
    
    int reachable_targets = 1;
    
    while (!q.empty()) {
        VertexId u = q.front();
        q.pop();
        
        for (size_t i = kao_[u]; i < kao_[u + 1]; ++i) {
            VertexId v = fo_[i];
            if (spot[v] == 0) {
                spot[v] = 1;
                q.push(v);
                if (target_vertices_[v] == 1) reachable_targets++;
            }
        }
    }
    
    return reachable_targets == target_count;
}

bool ReliabilityGraph::checkDistanceConstraint(int max_distance) const {
    size_t n = numVertices();
    int inf = max_distance + 100;
    std::vector<std::vector<int>> dist(n, std::vector<int>(n, inf));
    
    for (size_t i = 0; i < n; ++i) dist[i][i] = 0;
    
    for (size_t u = 0; u < n; ++u) {
        for (size_t k = kao_[u]; k < kao_[u + 1]; ++k) {
            VertexId v = fo_[k];
            dist[u][v] = 1;
        }
    }
    
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (dist[i][k] < inf && dist[k][j] < inf) {
                    dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }
    
    for (size_t i = 0; i < n; ++i) {
        if (target_vertices_[i]) {
            for (size_t j = i + 1; j < n; ++j) {
                if (target_vertices_[j] && dist[i][j] > max_distance) {
                    return false;
                }
            }
        }
    }
    
    return true;
}

// -----------------------------------------------------------------------------
// Core Reliability Calculation
// -----------------------------------------------------------------------------

ReliabilityResult ReliabilityGraph::calculateReliabilityWithDiameter(int diameter) const {
    long long recursions = 0;
    double total_rel = 0.0;
    
    std::function<void(ReliabilityGraph, int, double)> solve = 
        [&](ReliabilityGraph g, int variant, double current_rel) {
        recursions++;
        
        if (!variant && !g.isKConnected()) return;
        if (variant && !g.checkDistanceConstraint(diameter)) return;
        
        auto opt_edge = g.findLastUnreliableEdge();
        if (!opt_edge) {
            total_rel += current_rel;
            return;
        }
        
        EdgeId edge_idx = *opt_edge;
        double p = g.p_array_[edge_idx];
        g.p_array_[edge_idx] = 1.0; 
        
        auto [from, reverse_idx] = g.findReverseEdgeIndices(edge_idx);
        if (reverse_idx != -1) g.p_array_[reverse_idx] = 1.0;
        
        solve(g, 1, current_rel * p);
        
        if (reverse_idx != -1) g.p_array_[reverse_idx] = p;
        g.p_array_[edge_idx] = p;
        
        ReliabilityGraph t = g.removeEdge(g.fo_[edge_idx], from);
        solve(t, 0, current_rel * (1.0 - p));
    };
    
    solve(*this, 0, 1.0);
    
    return ReliabilityResult(total_rel, recursions, 0.0);
}

ReliabilityResult ReliabilityGraph::calculateReliabilityBetweenVertices(VertexId source, VertexId target, int diameter) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    LOG_DEBUG("Starting calculateReliabilityBetweenVertices: s={}, t={}, diameter={}", source, target, diameter);
    
    long long recursions = 0;
    double total_rel = 0.0;
    
    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    static thread_local long long last_logged_recursions = 0;
    constexpr long long LOG_INTERVAL = 100000; // Log every 100k recursions
    #endif
    
    std::function<void(ReliabilityGraph, int, double)> solve = 
        [&](ReliabilityGraph g, int variant, double current_rel) {
        recursions++;
        
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        if (recursions - last_logged_recursions >= LOG_INTERVAL) {
            LOG_DEBUG("Recursion progress: {} recursions, current_reliability={}", recursions, total_rel);
            last_logged_recursions = recursions;
        }
        #endif
        
        if (!variant) {
            if (g.calculateDistance(source, target, diameter) > diameter) return;
        }
        
        auto opt_edge = g.findLastUnreliableEdge();
        if (!opt_edge) {
            total_rel += current_rel;
            return;
        }
        
        EdgeId edge_idx = *opt_edge;
        double p = g.p_array_[edge_idx];
        g.p_array_[edge_idx] = 1.0;
        
        auto [from, reverse_idx] = g.findReverseEdgeIndices(edge_idx);
        if (reverse_idx != -1) g.p_array_[reverse_idx] = 1.0;
        
        solve(g, 1, current_rel * p);
        
        if (reverse_idx != -1) g.p_array_[reverse_idx] = p;
        g.p_array_[edge_idx] = p;
        
        ReliabilityGraph t = g.removeEdge(g.fo_[edge_idx], from);
        solve(t, 0, current_rel * (1.0 - p));
    };
    
    solve(*this, 0, 1.0);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration<double>(end_time - start_time).count();
    
    LOG_INFO("calculateReliabilityBetweenVertices completed: reliability={}, recursions={}, time={}s",
             total_rel, recursions, execution_time);
    
    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    last_logged_recursions = 0; // Reset for next call
    #endif
    
    return ReliabilityResult(total_rel, recursions, execution_time);
}

// Cancela-Petingi path-based factoring (Method 4)
// Optimized with global ISPT (compares feasible paths only) + per-pivot ISPT
ReliabilityResult ReliabilityGraph::calculateReliabilityCancelaPetingi(VertexId source,
                                                                       VertexId target,
                                                                       int diameter) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    long long recursions = 0;

    auto paths = PathEnumerator::enumeratePaths(*this, source, target, diameter);
    if (paths.empty()) {
        auto end_time = std::chrono::high_resolution_clock::now();
        double execution_time = std::chrono::duration<double>(end_time - start_time).count();
        return ReliabilityResult(0.0, 0, execution_time);
    }

    CancelaPetingiState state(paths, *this);
    
    state.applyGlobalISPT();

    std::function<double(CancelaPetingiState&)> facto = 
        [&](CancelaPetingiState& st) -> double {
        
        if (st.isTerminalSuccess()) return 1.0;
        if (st.isTerminalFail()) return 0.0;
        
        EdgeId e = st.selectPivotEdge();
        if (e == static_cast<EdgeId>(-1)) {
            return st.isTerminalSuccess() ? 1.0 : 0.0;
        }
        
        recursions++;
        
        st.applyISPT(e);
        double re = st.edge_reliability_[e];
        
        CancelaPetingiState st_contract(st);
        st_contract.applyContract(e);
        st_contract.edge_reliability_[e] = 1.0;
        double RContract = facto(st_contract);
        
        CancelaPetingiState st_delete(st);
        st_delete.applyDelete(e);
        st_delete.edge_reliability_.erase(e);
        st_delete.P_of_edge_.erase(e);
        double RDelete = facto(st_delete);
        
        return re * RContract + (1.0 - re) * RDelete;
    };

    double reliability = facto(state);

    auto end_time = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration<double>(end_time - start_time).count();

    LOG_INFO("calculateReliabilityCancelaPetingi completed: reliability={}, recursions={}, time={}s",
             reliability, recursions, execution_time);

    return ReliabilityResult(reliability, recursions, execution_time);
}

// Multi-diameter Cancela-Petingi: computes reliability CDF for multiple diameters in ONE pass
// Returns cumulative reliabilities: results[i] = P(path exists with length <= lowerBound + i)
std::pair<std::vector<double>, long long> ReliabilityGraph::calculateReliabilityCancelaPetingiMulti(
    VertexId source, VertexId target, int lowerBound, int upperBound) const {
    
    int range = upperBound - lowerBound + 1;
    std::vector<double> results(range, 0.0);
    long long recursions = 0;

    PROFILE_START(path_enum);
    auto paths = PathEnumerator::enumeratePaths(*this, source, target, upperBound);
    PROFILE_END(path_enum, path_enumeration_time);
    
    if (paths.empty()) {
        return {results, 0};
    }
    PROFILE_ADD_PATHS(paths.size());

    CancelaPetingiStateMulti state(paths, *this, upperBound);
    
    PROFILE_START(global_ispt);
    state.applyGlobalISPT();
    PROFILE_END(global_ispt, global_ispt_time);

    PROFILE_START(factoring);
    // factoMulti returns CDF: cdf[d] = P(exists path of length <= d)
    // NO early termination on "success" - we must compute CDF for all lengths
    std::function<std::vector<double>(CancelaPetingiStateMulti&)> factoMulti = 
        [&](CancelaPetingiStateMulti& st) -> std::vector<double> {
        
        std::vector<double> cdf(upperBound + 1, 0.0);
        
        // Terminal case: no more edges to factor
        // Set CDF based on current perfect paths
        EdgeId e = st.selectPivotEdge();
        if (e == static_cast<EdgeId>(-1)) {
            // No unreliable edges left - check which lengths have perfect paths
            for (int d = 1; d <= upperBound; ++d) {
                if (st.isTerminalSuccess(d)) {
                    cdf[d] = 1.0;
                }
            }
            // Ensure cumulative property
            for (int d = 1; d <= upperBound; ++d) {
                if (cdf[d-1] > cdf[d]) cdf[d] = cdf[d-1];
            }
            return cdf;
        }
        
        // Terminal case: no feasible paths
        if (st.allPathsInfeasible()) {
            return cdf; // all zeros
        }
        
        recursions++;
        
        st.applyISPT(e);
        double re = st.edge_reliability_[e];
        
        CancelaPetingiStateMulti st_contract(st);
        st_contract.applyContract(e);
        st_contract.edge_reliability_[e] = 1.0;
        auto cdf_contract = factoMulti(st_contract);
        
        CancelaPetingiStateMulti st_delete(st);
        st_delete.applyDelete(e);
        st_delete.edge_reliability_.erase(e);
        st_delete.P_of_edge_.erase(e);
        auto cdf_delete = factoMulti(st_delete);
        
        for (int d = 0; d <= upperBound; ++d) {
            cdf[d] = re * cdf_contract[d] + (1.0 - re) * cdf_delete[d];
        }
        return cdf;
    };

    auto cdf = factoMulti(state);
    PROFILE_END(factoring, factoring_time);
    
    // Extract results for requested range [lowerBound, upperBound]
    for (int i = 0; i < range; ++i) {
        int d = lowerBound + i;
        results[i] = (d <= upperBound) ? cdf[d] : cdf[upperBound];
    }

    return {results, recursions};
}

// Modified factoring - computes reliabilities for ALL diameters from lowerBound to upperBound in ONE pass
// This is the key optimization of M-Decomposition (Factoring2VertM from legacy code)
std::pair<std::vector<double>, long long> ReliabilityGraph::calculateReliabilityMultipleDiameters(
    VertexId source, VertexId target, int lowerBound, int upperBound) const {
    
    int range = upperBound - lowerBound + 1;
    std::vector<double> results(range, 0.0);
    long long recursions = 0;
    
    // Modified factoring: if distance > d but <= upperBound, continue with d+1
    std::function<void(ReliabilityGraph, int, int, double)> solveM = 
        [&](ReliabilityGraph g, int variant, int d, double current_rel) {
        recursions++;
        
        int dist = g.calculateDistance(source, target, upperBound + 1);
        
        if (!variant && dist > d) {
            if (dist <= upperBound) {
                recursions--; // Don't count this as a recursion (matches legacy behavior)
                solveM(g, 0, d + 1, current_rel);
            }
            return; // distance > upperBound, add 0
        }
        
        auto opt_edge = g.findLastUnreliableEdge();
        if (!opt_edge) {
            // Add reliability to result for diameter d
            if (d >= lowerBound && d <= upperBound) {
                results[d - lowerBound] += current_rel;
            }
            return;
        }
        
        EdgeId edge_idx = *opt_edge;
        double p = g.p_array_[edge_idx];
        g.p_array_[edge_idx] = 1.0;
        
        auto [from, reverse_idx] = g.findReverseEdgeIndices(edge_idx);
        if (reverse_idx != -1) g.p_array_[reverse_idx] = 1.0;
        
        solveM(g, 1, d, current_rel * p);
        
        if (reverse_idx != -1) g.p_array_[reverse_idx] = p;
        g.p_array_[edge_idx] = p;
        
        ReliabilityGraph t = g.removeEdge(g.fo_[edge_idx], from);
        solveM(t, 0, d, current_rel * (1.0 - p));
    };
    
    solveM(*this, 0, lowerBound, 1.0);
    
    // Convert to cumulative (results[i] = P(distance <= lowerBound + i))
    for (int i = 1; i < range; ++i) {
        results[i] += results[i - 1];
    }
    
    return {results, recursions};
}

// -----------------------------------------------------------------------------
// Decomposition Methods
// -----------------------------------------------------------------------------

// Helper for DFS
void searchVertices(const Graph& g, 
                   Graph::VertexId v, 
                   Graph::VertexId p, // Parent
                   int& l, // Counter
                   int& r, // Block counter
                   std::vector<int>& df_number, 
                   std::vector<int>& low, 
                   std::vector<int>& tree_edges, 
                   std::vector<int>& stack, 
                   std::vector<int>& result) {
    int last = 0;
    df_number[v] = l;
    low[v] = l++;

    for (size_t i = g.kao_[v]; i < g.kao_[v + 1]; ++i) {
        Graph::VertexId u = g.fo_[i];
        if (u == p) continue;

        // If not visited or visited but not tree edge (back edge)
        if (!df_number[u] || (df_number[u] < df_number[v] && !tree_edges[i])) {
            last = stack.size();
            stack.push_back(i);
            // Find reverse edge index
            auto rev_idx = g.findEdge(u, v);
            if (rev_idx != -1) stack.push_back(rev_idx);
        }

        if (!df_number[u]) {
            tree_edges[i] = 1;
            auto rev_idx = g.findEdge(u, v);
            if (rev_idx != -1) tree_edges[rev_idx] = 1;

            searchVertices(g, u, v, l, r, df_number, low, tree_edges, stack, result);

            if (low[u] >= df_number[v]) {
                for (size_t j = last; j < stack.size(); ++j) {
                    result[stack[j]] = r;
                }
                r++;
                stack.resize(last);
            }
            low[v] = std::min(low[v], low[u]);
        } else if (!tree_edges[i] && df_number[u] < df_number[v]) {
            low[v] = std::min(low[v], df_number[u]);
        }
    }
}

std::vector<int> ReliabilityGraph::decomposeIntoBlocks() const {
    size_t n = numVertices();
    size_t m = numEdges();
    std::vector<int> df_number(n, 0);
    std::vector<int> low(n, 0);
    std::vector<int> tree_edges(m, 0);
    std::vector<int> stack;
    int r = 1; 
    int l = 1;
    std::vector<int> result(m + 1, 0);

    if (n > 0) {
        // Start DFS from first vertex (0)
        searchVertices(*this, 0, -1, l, r, df_number, low, tree_edges, stack, result);
    }
    
    result[m] = r - 1; // Store number of blocks
    return result;
}

bool ReliabilityGraph::checkConnectivityWithoutBlock(int block_number, const std::vector<int>& block_decomposition) const {
    // Logic from ConnectivityWithoutBlock
    size_t n = numVertices();
    if (n == 0) return true;
    
    // Find first target vertex
    VertexId start_node = -1;
    int target_count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (target_vertices_[i]) {
            target_count++;
            if (start_node == -1) start_node = i;
        }
    }
    
    if (start_node == -1) return true; // No targets
    
    std::vector<int> spot(n, 0);
    std::vector<VertexId> q;
    q.push_back(start_node);
    spot[start_node] = 1;
    
    int reachable_targets = 1;
    
    size_t head = 0;
    while(head < q.size()){
        VertexId u = q[head++];
        for (size_t i = kao_[u]; i < kao_[u+1]; ++i) {
            VertexId v = fo_[i];
            // If edge i is NOT in block_number
            if (!spot[v] && block_decomposition[i] != block_number) {
                spot[v] = 1;
                q.push_back(v);
                if (target_vertices_[v]) reachable_targets++;
            }
        }
    }
    
    return reachable_targets == target_count;
}

std::vector<int> ReliabilityGraph::decomposeIntoKBlocks() const {
    std::vector<int> s = decomposeIntoBlocks();
    int num_blocks = s.back();
    if (num_blocks <= 1) return s; // 1 block or 0

    std::vector<int> target_blocks(num_blocks + 1, 0);
    size_t n = numVertices();
    
    // Identify blocks containing target vertices
    for (size_t i = 0; i < n; ++i) {
        if (target_vertices_[i]) {
            bool same_block = true;
            int k = -1;
            
            // Check edges of i
            for (size_t j = kao_[i]; j < kao_[i+1]; ++j) {
                if (k == -1) k = s[j];
                else if (s[j] != k) same_block = false;
            }
            
            if (same_block && k != -1) {
                target_blocks[k] = 1;
            }
        }
    }
    
    // Check connectivity without each block
    for (int i = 1; i <= num_blocks; ++i) {
        if (target_blocks[i] == 0 && !checkConnectivityWithoutBlock(i, s)) {
            target_blocks[i] = 1;
        }
    }
    
    // Refine targets? (Legacy logic: if vertex is not target but crucial for block connection?)
    // Legacy: "TargetOrder[!y!] = *amount of targetBlocks*;" - not clear but code follows.
    // Legacy code modifies this->Targets (mutable!) inside DecomposeOnBlocksK.
    // We cannot modify const member. But we can work with local copy if needed.
    // However, target_vertices_ is not mutable here. 
    // Legacy comment says "TargetOrder...". The loop actually checks if a non-target vertex 
    // connects only target blocks, then marks it as target?
    // For now, we skip modifying targets as it might be a side effect we don't want in a const method.
    
    // Remap block numbers
    std::vector<int> new_numbers(num_blocks + 1);
    for (int i = 0; i <= num_blocks; ++i) new_numbers[i] = i;
    
    int k = 0;
    for (int i = 1; i <= num_blocks; ++i) {
        if (target_blocks[i] == 0) {
            k++;
            new_numbers[i] = 0;
            for (int j = i + 1; j <= num_blocks; ++j) {
                new_numbers[j]--;
            }
        }
    }
    
    for (size_t i = 0; i < s.size() - 1; ++i) {
        s[i] = new_numbers[s[i]];
    }
    s.back() -= k;
    
    return s;
}

// -----------------------------------------------------------------------------
// M-Decomposition Support
// -----------------------------------------------------------------------------

struct BlockGraphNode {
    int original_block_id;
    std::vector<int> vertices_original_ids;
};

struct BlockGraphEdge {
    int from_block_node_idx;
    int to_block_node_idx;
    int articulation_point_original_id;
};

// Helper to find path in block graph
std::vector<int> findPathInBlockGraph(
    int start_block_node_idx,
    int end_block_node_idx,
    int num_block_nodes,
    const std::vector<BlockGraphEdge>& block_edges,
    std::vector<int>& out_path_articulation_points
) {
    if (start_block_node_idx == end_block_node_idx) {
        out_path_articulation_points.clear();
        return {start_block_node_idx};
    }

    std::vector<std::vector<std::pair<int, int>>> adj(num_block_nodes);
    for (const auto& edge : block_edges) {
        adj[edge.from_block_node_idx].push_back({edge.to_block_node_idx, edge.articulation_point_original_id});
        adj[edge.to_block_node_idx].push_back({edge.from_block_node_idx, edge.articulation_point_original_id});
    }

    std::queue<std::pair<int, std::vector<int>>> q; // node, path
    std::queue<std::vector<int>> q_aps; // articulation points path
    
    q.push({start_block_node_idx, {start_block_node_idx}});
    q_aps.push({});
    
    std::vector<bool> visited(num_block_nodes, false);
    visited[start_block_node_idx] = true;
    
    while(!q.empty()) {
        auto [u, path] = q.front(); q.pop();
        auto aps = q_aps.front(); q_aps.pop();
        
        if (u == end_block_node_idx) {
            out_path_articulation_points = aps;
            return path;
        }
        
        for (auto& edge : adj[u]) {
            int v = edge.first;
            int ap = edge.second;
            
            if (!visited[v]) {
                visited[v] = true;
                
                std::vector<int> new_path = path;
                new_path.push_back(v);
                q.push({v, new_path});
                
                std::vector<int> new_aps = aps;
                new_aps.push_back(ap);
                q_aps.push(new_aps);
            }
        }
    }
    
    return {};
}

std::vector<int> ReliabilityGraph::getBlocksContainingVertex(VertexId v, const std::vector<int>& decomposition) const {
    std::vector<int> blocks;
    if (v < 0 || v >= static_cast<VertexId>(numVertices())) return blocks;
    
    int num_blocks = decomposition.back();
    std::vector<bool> seen(num_blocks + 1, false);
    
    for (size_t i = kao_[v]; i < kao_[v+1]; ++i) {
        if (i < decomposition.size() - 1) {
            int b = decomposition[i];
            if (b > 0 && !seen[b]) {
                seen[b] = true;
                blocks.push_back(b);
            }
        }
    }
    return blocks;
}

int ReliabilityGraph::calculateMinBlockDiameter(
    int block_id,
    const std::vector<int>& decomposition,
    VertexId entry,
    VertexId exit
) const {
    if (entry == exit) return 0;
    
    ReliabilityGraph block_graph;
    std::vector<int> map_new_to_orig, map_orig_to_new;
    getBlockGraphAndMap(block_id, decomposition, block_graph, map_new_to_orig, map_orig_to_new);
    
    if (entry >= static_cast<VertexId>(map_orig_to_new.size()) ||
        exit >= static_cast<VertexId>(map_orig_to_new.size())) {
        return -1;
    }
    
    int s_local = map_orig_to_new[entry];
    int t_local = map_orig_to_new[exit];
    
    if (s_local == -1 || t_local == -1) {
        return -1;
    }
    
    // BFS to find shortest path
    size_t n = block_graph.numVertices();
    std::vector<int> dist(n, -1);
    std::queue<int> q;
    q.push(s_local);
    dist[s_local] = 0;
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        if (u == t_local) {
            return dist[u];
        }
        
        for (size_t i = block_graph.kao_[u]; i < block_graph.kao_[u + 1]; ++i) {
            int v = block_graph.fo_[i];
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    
    return -1;  // unreachable
}

void ReliabilityGraph::getBlockGraphAndMap(
    int block_id,
    const std::vector<int>& decomposition,
    ReliabilityGraph& out_graph,
    std::vector<int>& map_new_to_orig,
    std::vector<int>& map_orig_to_new
) const {
    size_t n = numVertices();
    map_orig_to_new.assign(n, -1);
    map_new_to_orig.clear();
    
    std::set<VertexId> block_vertices;
    std::vector<std::pair<VertexId, VertexId>> edges;
    std::vector<Probability> probs;
    
    // Identify vertices and edges in this block
    for (size_t v = 0; v < n; ++v) {
        for (size_t i = kao_[v]; i < kao_[v+1]; ++i) {
            if (i < decomposition.size() - 1 && decomposition[i] == block_id) {
                VertexId u = fo_[i];
                if (v < u) { // Add each edge once
                    block_vertices.insert(v);
                    block_vertices.insert(u);
                    edges.push_back({v, u});
                    probs.push_back(p_array_[i]);
                }
            }
        }
    }
    
    // Build mappings
    int new_id = 0;
    for (VertexId v : block_vertices) {
        map_orig_to_new[v] = new_id;
        map_new_to_orig.push_back(v);
        new_id++;
    }
    
    // Build CSR
    int num_new_v = new_id;
    if (num_new_v == 0) {
        out_graph = ReliabilityGraph({0}, {}, {}, {});
        return;
    }
    
    std::vector<std::vector<int>> adj(num_new_v);
    std::vector<Probability> p_array_new; // This assumes order is preserved or reconstructed?
    // Reconstructing CSR is tricky with edge probabilities if not sequential.
    // Let's just build adj list first, then flatten.
    
    struct EdgeInfo { int to; double p; };
    std::vector<std::vector<EdgeInfo>> adj_info(num_new_v);
    
    for (size_t k = 0; k < edges.size(); ++k) {
        int u_orig = edges[k].first;
        int v_orig = edges[k].second;
        double p = probs[k];
        
        int u_new = map_orig_to_new[u_orig];
        int v_new = map_orig_to_new[v_orig];
        
        adj_info[u_new].push_back({v_new, p});
        adj_info[v_new].push_back({u_new, p});
    }
    
    std::vector<VertexId> new_kao(num_new_v + 1);
    new_kao[0] = 0;
    std::vector<VertexId> new_fo;
    std::vector<Probability> new_p_array;
    
    for (int i = 0; i < num_new_v; ++i) {
        // Sort neighbors to be deterministic?
        // Not strictly required for reliability but good for consistency.
        std::sort(adj_info[i].begin(), adj_info[i].end(), [](const EdgeInfo& a, const EdgeInfo& b){
            return a.to < b.to;
        });
        
        new_kao[i+1] = new_kao[i] + adj_info[i].size();
        for (const auto& e : adj_info[i]) {
            new_fo.push_back(e.to);
            new_p_array.push_back(e.p);
        }
    }
    
    std::vector<int> new_targets(num_new_v, 0);
    for (int i = 0; i < num_new_v; ++i) {
        int orig = map_new_to_orig[i];
        if (orig < target_vertices_.size()) {
            new_targets[i] = target_vertices_[orig];
        }
    }
    
    out_graph = ReliabilityGraph(std::move(new_kao), std::move(new_fo), std::move(new_p_array), std::move(new_targets));
}

// Recursive solver for block chain
std::vector<double> ReliabilityGraph::solveRecursiveForBlockChain(
    int current_idx,
    int entry_node_orig,
    int target_node_orig,
    int max_len,
    const std::vector<int>& decomposition,
    const std::vector<int>& block_ids,
    const std::vector<int>& aps_orig,
    long long& recursion_counter
) const {
    recursion_counter++;
    
    if (max_len < 0 || current_idx >= block_ids.size()) return {};
    
    int block_id = block_ids[current_idx];
    ReliabilityGraph block_graph;
    std::vector<int> map_new_to_orig, map_orig_to_new;
    
    getBlockGraphAndMap(block_id, decomposition, block_graph, map_new_to_orig, map_orig_to_new);
    
    if (block_graph.numVertices() <= 1 || map_orig_to_new[entry_node_orig] == -1) {
        return std::vector<double>(max_len + 1, 0.0);
    }
    
    int s_local = map_orig_to_new[entry_node_orig];
    bool is_last = (current_idx == block_ids.size() - 1);
    
    // Check if target is in this block
    bool target_in_block = false;
    auto target_blocks = getBlocksContainingVertex(target_node_orig, decomposition);
    for(int b : target_blocks) if(b == block_id) target_in_block = true;
    
    // For large blocks, fall back to CPFM to avoid exponential blowup in modified factoring
    const bool use_cpfm_fallback = (block_graph.numEdges() > 30);

    if (is_last && target_in_block) {
        if (map_orig_to_new[target_node_orig] == -1) return std::vector<double>(max_len + 1, 0.0);
        int t_local = map_orig_to_new[target_node_orig];

        std::vector<double> r_cumulative;
        long long recs = 0;
        if (use_cpfm_fallback) {
            auto [rc, r] = block_graph.calculateReliabilityCancelaPetingiMulti(s_local, t_local, 0, max_len);
            r_cumulative = std::move(rc); recs = r;
        } else {
            auto [rc, r] = block_graph.calculateReliabilityMultipleDiameters(s_local, t_local, 0, max_len);
            r_cumulative = std::move(rc); recs = r;
        }
        recursion_counter += recs;

        return r_cumulative;
    }

    if (is_last && !target_in_block) return std::vector<double>(max_len + 1, 0.0);

    if (current_idx + 1 >= aps_orig.size()) return std::vector<double>(max_len + 1, 0.0);
    int exit_node_orig = aps_orig[current_idx + 1];

    if (map_orig_to_new[exit_node_orig] == -1) return std::vector<double>(max_len + 1, 0.0);
    int x_local = map_orig_to_new[exit_node_orig];

    // Calculate S->X in this block — for large blocks use CPFM fallback
    std::vector<double> r_cumulative_sx;
    long long recs_sx = 0;
    if (use_cpfm_fallback) {
        auto [rc, r] = block_graph.calculateReliabilityCancelaPetingiMulti(s_local, x_local, 0, max_len);
        r_cumulative_sx = std::move(rc); recs_sx = r;
    } else {
        auto [rc, r] = block_graph.calculateReliabilityMultipleDiameters(s_local, x_local, 0, max_len);
        r_cumulative_sx = std::move(rc); recs_sx = r;
    }
    recursion_counter += recs_sx;
    
    // Compute PMF from CDF
    std::vector<double> r_pmf_sx(max_len + 1);
    r_pmf_sx[0] = r_cumulative_sx[0];
    for (int d = 1; d <= max_len; ++d) {
        r_pmf_sx[d] = r_cumulative_sx[d] - r_cumulative_sx[d-1];
    }
    
    // Recursive call
    std::vector<double> r_cumulative_next = solveRecursiveForBlockChain(
        current_idx + 1, exit_node_orig, target_node_orig, 
        max_len, decomposition, block_ids, aps_orig, recursion_counter
    );
    
    // Convolve
    std::vector<double> result(max_len + 1, 0.0);
    for (int d = 0; d <= max_len; ++d) {
        double sum = 0.0;
        for (int k = 0; k <= d; ++k) {
            if (r_pmf_sx[k] > 0) {
                int remaining = d - k;
                if (remaining < r_cumulative_next.size()) {
                    sum += r_pmf_sx[k] * r_cumulative_next[remaining];
                }
            }
        }
        result[d] = sum;
    }
    for (int d = 1; d <= max_len; ++d) {
        if (result[d] < result[d-1]) result[d] = result[d-1];
    }
    return result;
}

// =============================================================================
// LEVEL 1: TRUE NESTED RECURSION (The "Naive" Method)
// =============================================================================
// This implements the theoretically correct but INEFFICIENT nested recursion.
// For each path length L in block_i, we make a RECURSIVE CALL to block_{i+1}.
// This causes redundant computation (the "memory effect") which is intentional
// for academic comparison purposes.
// =============================================================================

double ReliabilityGraph::solveNestedRecursive(
    int current_idx,
    int entry_node_orig,
    int target_node_orig,
    int remaining_diameter,
    const std::vector<int>& decomposition,
    const std::vector<int>& block_ids,
    const std::vector<int>& aps_orig,
    long long& recursion_counter
) const {
    recursion_counter++;
    
    // Base case: invalid parameters
    if (remaining_diameter < 0 || current_idx >= static_cast<int>(block_ids.size())) {
        return 0.0;
    }
    
    int block_id = block_ids[current_idx];
    ReliabilityGraph block_graph;
    std::vector<int> map_new_to_orig, map_orig_to_new;
    
    getBlockGraphAndMap(block_id, decomposition, block_graph, map_new_to_orig, map_orig_to_new);
    
    if (block_graph.numVertices() <= 1 || map_orig_to_new[entry_node_orig] == -1) {
        return 0.0;
    }
    
    int s_local = map_orig_to_new[entry_node_orig];
    bool is_last = (current_idx == static_cast<int>(block_ids.size()) - 1);
    
    // Check if target is in this block
    bool target_in_block = false;
    auto target_blocks = getBlocksContainingVertex(target_node_orig, decomposition);
    for (int b : target_blocks) {
        if (b == block_id) target_in_block = true;
    }
    
    // LAST BLOCK: Return cumulative reliability P(distance <= remaining_diameter)
    if (is_last && target_in_block) {
        if (map_orig_to_new[target_node_orig] == -1) return 0.0;
        int t_local = map_orig_to_new[target_node_orig];
        
        // Standard factoring for single diameter value
        auto res = block_graph.calculateReliabilityBetweenVertices(s_local, t_local, remaining_diameter);
        recursion_counter += res.recursions;
        return res.reliability;
    }
    
    if (is_last && !target_in_block) return 0.0;
    
    // Not the last block: need to recurse to next block
    if (current_idx + 1 >= static_cast<int>(aps_orig.size())) return 0.0;
    int exit_node_orig = aps_orig[current_idx + 1];
    
    if (map_orig_to_new[exit_node_orig] == -1) return 0.0;
    int x_local = map_orig_to_new[exit_node_orig];
    
    // ==========================================================================
    // NESTED RECURSION: The "Memory Effect"
    // For each possible path length L in this block, we:
    //   1. Compute P(distance == L) = P(distance <= L) - P(distance <= L-1)
    //   2. Make a RECURSIVE CALL for the remaining chain with (remaining_diameter - L)
    //   3. Accumulate: total += P(exact L) * RecursiveResult
    // ==========================================================================
    
    double total_reliability = 0.0;
    
    // Iterate through all possible path lengths in current block
    for (int len = 0; len <= remaining_diameter; ++len) {
        // Calculate P(distance <= len) using standard factoring
        auto res_len = block_graph.calculateReliabilityBetweenVertices(s_local, x_local, len);
        recursion_counter += res_len.recursions;
        double val_len = res_len.reliability;
        
        // Calculate P(distance <= len-1) for computing exact probability
        double val_prev = 0.0;
        if (len > 0) {
            auto res_prev = block_graph.calculateReliabilityBetweenVertices(s_local, x_local, len - 1);
            recursion_counter += res_prev.recursions;
            val_prev = res_prev.reliability;
        }
        
        // P(distance == exactly len) = P(distance <= len) - P(distance <= len-1)
        double prob_exact = val_len - val_prev;
        
        if (prob_exact > 1e-15) {
            // CRITICAL: Recursive call happens HERE, INSIDE the loop
            // This triggers re-calculation of the tail for every length L
            // This is the "memory effect" - intentionally inefficient!
            double tail_reliability = solveNestedRecursive(
                current_idx + 1,
                exit_node_orig,
                target_node_orig,
                remaining_diameter - len,
                decomposition,
                block_ids,
                aps_orig,
                recursion_counter
            );
            
            total_reliability += prob_exact * tail_reliability;
        }
    }
    
    return total_reliability;
}

// -----------------------------------------------------------------------------
// Public Calculation Methods
// -----------------------------------------------------------------------------

ReliabilityResult ReliabilityGraph::calculateReliabilityWithMDecomposition(VertexId s, VertexId t, int d) const {
    LOG_INFO("Starting M-Decomposition: s={}, t={}, diameter={}", s, t, d);

    clock_t start = clock();
    long long recursion_counter = 0;

    LOG_DEBUG("Decomposing graph into blocks");
    std::vector<int> decomposition = decomposeIntoBlocks();
    int num_blocks = decomposition.back();
    LOG_DEBUG("Graph decomposed into {} blocks", num_blocks);
    
    // If single block, fall back to standard factoring
    if (num_blocks <= 1) {
        LOG_DEBUG("Single block detected, using standard factoring");
        return calculateReliabilityBetweenVertices(s, t, d);
    }
    
    // Build Block Graph
    std::vector<BlockGraphNode> block_nodes(num_blocks);
    std::map<int, int> block_id_to_idx;
    for(int i=0; i<num_blocks; ++i) {
        block_nodes[i].original_block_id = i + 1;
        block_id_to_idx[i+1] = i;
    }
    
    std::vector<BlockGraphEdge> block_edges;
    std::vector<bool> ap_processed(numVertices(), false);
    
    for (size_t v = 0; v < numVertices(); ++v) {
        auto blocks = getBlocksContainingVertex(v, decomposition);
        if (blocks.size() > 1) {
            if (ap_processed[v]) continue;
            ap_processed[v] = true;
            for(size_t i=0; i<blocks.size(); ++i) {
                for(size_t j=i+1; j<blocks.size(); ++j) {
                    if (block_id_to_idx.count(blocks[i]) && block_id_to_idx.count(blocks[j])) {
                        block_edges.push_back({
                            block_id_to_idx[blocks[i]],
                            block_id_to_idx[blocks[j]],
                            static_cast<int>(v)
                        });
                    }
                }
            }
        }
    }
    
    auto s_blocks = getBlocksContainingVertex(s, decomposition);
    auto t_blocks = getBlocksContainingVertex(t, decomposition);
    
    if (s_blocks.empty() || t_blocks.empty()) {
         double time = (double)(clock() - start) / CLOCKS_PER_SEC;
         return ReliabilityResult(0.0, 0, time);
    }
    
    int start_node_idx = block_id_to_idx[s_blocks[0]];
    int end_node_idx = block_id_to_idx[t_blocks[0]];
    
    std::vector<int> path_aps;
    std::vector<int> block_path_indices = findPathInBlockGraph(
        start_node_idx, end_node_idx, num_blocks, block_edges, path_aps
    );
    
    if (block_path_indices.empty()) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, 0, time);
    }
    
    std::vector<int> ordered_block_ids;
    for (int idx : block_path_indices) {
        ordered_block_ids.push_back(block_nodes[idx].original_block_id);
    }
    
    std::vector<int> final_aps(ordered_block_ids.size() + 1, 0);
    for(size_t i=0; i<path_aps.size(); ++i) final_aps[i+1] = path_aps[i];
    
    LOG_DEBUG("Solving block chain with {} blocks", ordered_block_ids.size());
    std::vector<double> results = solveRecursiveForBlockChain(
        0, s, t, d, decomposition, ordered_block_ids, final_aps, recursion_counter
    );
    
    double rel = (results.empty()) ? 0.0 : results[std::min((int)results.size()-1, d)];
    double time = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    LOG_INFO("M-Decomposition completed: reliability={}, recursions={}, time={}s", rel, recursion_counter, time);
    
    return ReliabilityResult(rel, recursion_counter, time);
}

// Iterative block chain solver for Method 5 (replaces recursive version).
// Processes blocks from last to first (backward pass), accumulating CDF via convolution.
// Advantages over recursive version:
//   - No stack frame overhead or vector-by-value returns
//   - Precomputed block graphs (each constructed once)
//   - Suffix sums of d_min computed upfront
//
// Parameters:
//   block_graphs   — preextracted subgraph for each block in chain
//   s_local[i]     — entry vertex in block i (local numbering)
//   x_local[i]     — exit vertex in block i  (local numbering; for last block = t_local)
//   d_min_per_block — shortest s_local→x_local distance in block i
//   gap            — D − sum(d_min): extra "budget" shared across all blocks
static std::vector<double> solveBlockChainIterative(
    const std::vector<ReliabilityGraph>& block_graphs,
    const std::vector<int>& s_local,
    const std::vector<int>& x_local,
    const std::vector<int>& d_min_per_block,
    int gap,
    long long& recursion_counter
) {
    int n = static_cast<int>(block_graphs.size());
    if (n == 0) return {};

    // Precompute suffix sums: suffix[i] = d_min[i] + ... + d_min[n-1]
    std::vector<int> suffix(n + 1, 0);
    for (int i = n - 1; i >= 0; --i)
        suffix[i] = suffix[i + 1] + d_min_per_block[i];

    std::vector<double> acc_cdf;  // CDF of blocks [i .. n-1], grows leftward

    for (int i = n - 1; i >= 0; --i) {
        int d_min_k    = d_min_per_block[i];
        int max_len_k  = d_min_k + gap;           // max this block can consume
        int max_total  = suffix[i] + gap;          // max total for blocks [i..n-1]

        // CPFM: CDF for path length in [d_min_k, max_len_k] within this block
        auto [r_partial, recs] = block_graphs[i].calculateReliabilityCancelaPetingiMulti(
            s_local[i], x_local[i], d_min_k, max_len_k);
        recursion_counter += recs;

        if (i == n - 1) {
            // Last block — initialise acc_cdf directly
            acc_cdf.assign(max_total + 1, 0.0);
            for (int d = d_min_k; d <= max_len_k && d <= max_total; ++d)
                acc_cdf[d] = r_partial[d - d_min_k];
        } else {
            // Convert partial CDF → PMF for this block
            std::vector<double> pmf(max_len_k + 1, 0.0);
            pmf[d_min_k] = r_partial[0];
            for (int d = d_min_k + 1; d <= max_len_k; ++d)
                pmf[d] = r_partial[d - d_min_k] - r_partial[d - d_min_k - 1];

            // Convolve: new_cdf[d] = sum_k pmf[k] * acc_cdf[d - k]
            std::vector<double> new_cdf(max_total + 1, 0.0);
            for (int d = d_min_k; d <= max_total; ++d) {
                double sum = 0.0;
                int k_hi = std::min(d, max_len_k);
                for (int k = d_min_k; k <= k_hi; ++k) {
                    if (pmf[k] > 0.0) {
                        int rem = d - k;
                        if (rem < static_cast<int>(acc_cdf.size()))
                            sum += pmf[k] * acc_cdf[rem];
                    }
                }
                new_cdf[d] = sum;
            }
            acc_cdf = std::move(new_cdf);
        }

        // Enforce cumulative (non-decreasing) property
        for (int d = 1; d <= max_total; ++d)
            if (acc_cdf[d] < acc_cdf[d - 1]) acc_cdf[d] = acc_cdf[d - 1];
    }

    return acc_cdf;
}


ReliabilityResult ReliabilityGraph::calculateReliabilityWithMDecompositionCPFM(VertexId s, VertexId t, int d) const {
    LOG_INFO("Starting M-Decomposition + CPFM (Method 5): s={}, t={}, diameter={}", s, t, d);
    
    PROFILE_RESET();
    auto total_start = std::chrono::high_resolution_clock::now();
    long long recursion_counter = 0;
    
    PROFILE_START(decomp);
    LOG_DEBUG("Decomposing graph into blocks");
    std::vector<int> decomposition = decomposeIntoBlocks();
    int num_blocks = decomposition.back();
    PROFILE_END(decomp, decomposition_time);
    LOG_DEBUG("Graph decomposed into {} blocks", num_blocks);
    
#ifdef PROFILE_METHOD5
    g_profile.num_blocks = num_blocks;
#endif
    
    if (num_blocks <= 1) {
        LOG_DEBUG("Single block detected, using CPFM directly");
        return calculateReliabilityCancelaPetingi(s, t, d);
    }
    
    PROFILE_START(block_graph);
    std::vector<BlockGraphNode> block_nodes(num_blocks);
    std::map<int, int> block_id_to_idx;
    for (int i = 0; i < num_blocks; ++i) {
        block_nodes[i].original_block_id = i + 1;
        block_id_to_idx[i + 1] = i;
    }
    
    std::vector<BlockGraphEdge> block_edges;
    std::vector<bool> ap_processed(numVertices(), false);
    
    for (size_t v = 0; v < numVertices(); ++v) {
        auto blocks = getBlocksContainingVertex(v, decomposition);
        if (blocks.size() > 1) {
            if (ap_processed[v]) continue;
            ap_processed[v] = true;
            for (size_t i = 0; i < blocks.size(); ++i) {
                for (size_t j = i + 1; j < blocks.size(); ++j) {
                    if (block_id_to_idx.count(blocks[i]) && block_id_to_idx.count(blocks[j])) {
                        block_edges.push_back({
                            block_id_to_idx[blocks[i]],
                            block_id_to_idx[blocks[j]],
                            static_cast<int>(v)
                        });
                    }
                }
            }
        }
    }
    PROFILE_END(block_graph, block_graph_build_time);
    
    auto s_blocks = getBlocksContainingVertex(s, decomposition);
    auto t_blocks = getBlocksContainingVertex(t, decomposition);
    
    if (s_blocks.empty() || t_blocks.empty()) {
        auto total_end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(total_end - total_start).count();
        return ReliabilityResult(0.0, 0, time);
    }
    
    PROFILE_START(path_find);
    int start_node_idx = block_id_to_idx[s_blocks[0]];
    int end_node_idx = block_id_to_idx[t_blocks[0]];
    
    std::vector<int> path_aps;
    std::vector<int> block_path_indices = findPathInBlockGraph(
        start_node_idx, end_node_idx, num_blocks, block_edges, path_aps
    );
    PROFILE_END(path_find, path_finding_time);
    
    if (block_path_indices.empty()) {
        auto total_end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(total_end - total_start).count();
        return ReliabilityResult(0.0, 0, time);
    }
    
    std::vector<int> ordered_block_ids;
    for (int idx : block_path_indices) {
        ordered_block_ids.push_back(block_nodes[idx].original_block_id);
    }

    // If any block in the chain is too large for per-block CPFM, fall back to global CPFM (m4)
    // Threshold: >50 edges means per-block CPFM becomes exponentially expensive
    const int CPFM_BLOCK_EDGE_THRESHOLD = 50;
    for (int bid : ordered_block_ids) {
        ReliabilityGraph blk;
        std::vector<int> m2o, m2n;
        getBlockGraphAndMap(bid, decomposition, blk, m2o, m2n);
        if (static_cast<int>(blk.numEdges()) > CPFM_BLOCK_EDGE_THRESHOLD) {
            LOG_DEBUG("Block {} has {} edges > threshold, falling back to global CPFM (m4)",
                      bid, blk.numEdges());
            return calculateReliabilityCancelaPetingi(s, t, d);
        }
    }

    std::vector<int> final_aps(ordered_block_ids.size() + 1, 0);
    final_aps[0] = s;  // source vertex as first entry point
    for (size_t i = 0; i < path_aps.size(); ++i) {
        final_aps[i + 1] = path_aps[i];
    }

    // Calculate d_min for each block (formula 2.3 from thesis)
    std::vector<int> d_min_per_block(ordered_block_ids.size());
    int total_d_min = 0;
    
    for (size_t i = 0; i < ordered_block_ids.size(); ++i) {
        int block_id = ordered_block_ids[i];
        int entry = final_aps[i];
        int exit;
        
        // Determine exit vertex for this block
        if (i == ordered_block_ids.size() - 1) {
            // Last block: exit is target vertex
            exit = t;
        } else {
            exit = final_aps[i + 1];
        }
        
        int d_min_k = calculateMinBlockDiameter(block_id, decomposition, entry, exit);
        if (d_min_k < 0) {
            // Unreachable - should not happen in a valid block path
            auto total_end = std::chrono::high_resolution_clock::now();
            double time = std::chrono::duration<double>(total_end - total_start).count();
            return ReliabilityResult(0.0, 0, time);
        }
        d_min_per_block[i] = d_min_k;
        total_d_min += d_min_k;
    }
    
    int gap = d - total_d_min;
    LOG_DEBUG("Gap optimization: D={}, total_d_min={}, gap={}", d, total_d_min, gap);

    if (gap < 0) {
        auto total_end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(total_end - total_start).count();
        LOG_INFO("M-Decomposition + CPFM: gap<0, returning 0");
        return ReliabilityResult(0.0, 0, time);
    }

    // Precompute block subgraphs and local vertex indices once — avoids repeated
    // extraction inside the chain solver and eliminates temporary allocations per block.
    int chain_len = static_cast<int>(ordered_block_ids.size());
    std::vector<ReliabilityGraph> chain_block_graphs(chain_len);
    std::vector<int>              chain_s_local(chain_len);
    std::vector<int>              chain_x_local(chain_len);

    for (int i = 0; i < chain_len; ++i) {
        std::vector<int> m2o, o2n;
        getBlockGraphAndMap(ordered_block_ids[i], decomposition,
                            chain_block_graphs[i], m2o, o2n);

        int entry = final_aps[i];
        int exit_v = (i == chain_len - 1) ? t : final_aps[i + 1];

        if (entry  >= static_cast<int>(o2n.size()) || o2n[entry]  == -1 ||
            exit_v >= static_cast<int>(o2n.size()) || o2n[exit_v] == -1) {
            auto total_end = std::chrono::high_resolution_clock::now();
            double time = std::chrono::duration<double>(total_end - total_start).count();
            return ReliabilityResult(0.0, 0, time);
        }
        chain_s_local[i] = o2n[entry];
        chain_x_local[i] = o2n[exit_v];
    }

    LOG_DEBUG("Solving block chain ({} blocks, gap={}) using iterative CPFM solver",
              chain_len, gap);

    std::vector<double> results = solveBlockChainIterative(
        chain_block_graphs, chain_s_local, chain_x_local,
        d_min_per_block, gap, recursion_counter
    );
    
    // The result is indexed by total path length, so we need results[d]
    double rel = (results.empty() || d >= static_cast<int>(results.size())) ? 0.0 : results[d];
    auto total_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(total_end - total_start).count();

#ifdef PROFILE_METHOD5
    g_profile.total_time = time;
    g_profile.total_recursions = recursion_counter;
    PROFILE_PRINT();
#endif
    
    LOG_INFO("M-Decomposition + CPFM completed: reliability={}, recursions={}, time={}s", rel, recursion_counter, time);
    
    return ReliabilityResult(rel, recursion_counter, time);
}

ReliabilityResult ReliabilityGraph::calculateReliabilityWithRecursiveDecomposition(VertexId s, VertexId t, int d) const {
    // ==========================================================================
    // LEVEL 1: Recursive Decomposition (The "Naive" Method)
    // ==========================================================================
    // This method uses TRUE NESTED RECURSION where the solver for Block(i) makes
    // recursive calls to Block(i+1) INSIDE its loop over path lengths.
    // This is intentionally INEFFICIENT for academic comparison purposes.
    // ==========================================================================
    LOG_INFO("Starting Recursive Decomposition: s={}, t={}, diameter={}", s, t, d);
    clock_t start = clock();
    long long recursion_counter = 0;

    LOG_DEBUG("Decomposing graph into blocks");
    std::vector<int> decomposition = decomposeIntoBlocks();
    int num_blocks = decomposition.back();
    LOG_DEBUG("Graph decomposed into {} blocks", num_blocks);
    
    // If single block, fall back to standard factoring
    if (num_blocks <= 1) {
        LOG_DEBUG("Single block detected, using standard factoring");
        return calculateReliabilityBetweenVertices(s, t, d);
    }
    
    // Build Block Graph
    std::vector<BlockGraphNode> block_nodes(num_blocks);
    std::map<int, int> block_id_to_idx;
    for (int i = 0; i < num_blocks; ++i) {
        block_nodes[i].original_block_id = i + 1;
        block_id_to_idx[i + 1] = i;
    }
    
    std::vector<BlockGraphEdge> block_edges;
    std::vector<bool> ap_processed(numVertices(), false);
    
    for (size_t v = 0; v < numVertices(); ++v) {
        auto blocks = getBlocksContainingVertex(v, decomposition);
        if (blocks.size() > 1) {
            if (ap_processed[v]) continue;
            ap_processed[v] = true;
            for (size_t i = 0; i < blocks.size(); ++i) {
                for (size_t j = i + 1; j < blocks.size(); ++j) {
                    if (block_id_to_idx.count(blocks[i]) && block_id_to_idx.count(blocks[j])) {
                        block_edges.push_back({
                            block_id_to_idx[blocks[i]],
                            block_id_to_idx[blocks[j]],
                            static_cast<int>(v)
                        });
                    }
                }
            }
        }
    }
    
    auto s_blocks = getBlocksContainingVertex(s, decomposition);
    auto t_blocks = getBlocksContainingVertex(t, decomposition);
    
    if (s_blocks.empty() || t_blocks.empty()) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, 0, time);
    }
    
    int start_node_idx = block_id_to_idx[s_blocks[0]];
    int end_node_idx = block_id_to_idx[t_blocks[0]];
    
    std::vector<int> path_aps;
    std::vector<int> block_path_indices = findPathInBlockGraph(
        start_node_idx, end_node_idx, num_blocks, block_edges, path_aps
    );
    
    if (block_path_indices.empty()) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, 0, time);
    }
    
    std::vector<int> ordered_block_ids;
    for (int idx : block_path_indices) {
        ordered_block_ids.push_back(block_nodes[idx].original_block_id);
    }
    
    std::vector<int> final_aps(ordered_block_ids.size() + 1, 0);
    for (size_t i = 0; i < path_aps.size(); ++i) {
        final_aps[i + 1] = path_aps[i];
    }
    
    // Use TRUE NESTED RECURSION (not convolution!)
    // The recursive call happens INSIDE the loop over path lengths
    LOG_DEBUG("Solving nested recursive with {} blocks", ordered_block_ids.size());
    double reliability = solveNestedRecursive(
        0, s, t, d, decomposition, ordered_block_ids, final_aps, recursion_counter
    );
    
    double time = (double)(clock() - start) / CLOCKS_PER_SEC;
    LOG_INFO("Recursive Decomposition completed: reliability={}, recursions={}, time={}s", 
             reliability, recursion_counter, time);
    
    return ReliabilityResult(reliability, recursion_counter, time);
}

ReliabilityResult ReliabilityGraph::calculateReliabilityWithDecomposition(VertexId s, VertexId t, int UpperBound) const {
    // Uses block decomposition with simple (one-diameter-at-a-time) factoring per block.
    // Block chain is found via BFS on the block graph from s_block to t_block,
    // so this method correctly handles any (s,t) pair regardless of target_vertices_.
    LOG_INFO("Starting Simple Factoring with Decomposition: s={}, t={}, diameter={}", s, t, UpperBound);
    clock_t start = clock();
    long long recursion_counter = 0;

    LOG_DEBUG("Decomposing graph into blocks");
    std::vector<int> decomposition = decomposeIntoBlocks();
    int num_blocks = decomposition.back();
    LOG_DEBUG("Graph decomposed into {} blocks", num_blocks);

    // If single block, fall back to standard factoring
    if (num_blocks <= 1) {
        LOG_DEBUG("Single block detected, using standard factoring");
        return calculateReliabilityBetweenVertices(s, t, UpperBound);
    }

    // Build block graph (same as method 3)
    std::vector<BlockGraphNode> block_nodes(num_blocks);
    std::map<int, int> block_id_to_idx;
    for (int i = 0; i < num_blocks; ++i) {
        block_nodes[i].original_block_id = i + 1;
        block_id_to_idx[i + 1] = i;
    }

    std::vector<BlockGraphEdge> block_edges;
    std::vector<bool> ap_processed(numVertices(), false);
    for (size_t v = 0; v < numVertices(); ++v) {
        auto blocks = getBlocksContainingVertex(v, decomposition);
        if (blocks.size() > 1) {
            if (ap_processed[v]) continue;
            ap_processed[v] = true;
            for (size_t i = 0; i < blocks.size(); ++i)
                for (size_t j = i + 1; j < blocks.size(); ++j)
                    if (block_id_to_idx.count(blocks[i]) && block_id_to_idx.count(blocks[j]))
                        block_edges.push_back({block_id_to_idx[blocks[i]],
                                               block_id_to_idx[blocks[j]],
                                               static_cast<int>(v)});
        }
    }

    auto s_blocks = getBlocksContainingVertex(s, decomposition);
    auto t_blocks = getBlocksContainingVertex(t, decomposition);
    if (s_blocks.empty() || t_blocks.empty()) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, 0, time);
    }

    int start_block_idx = block_id_to_idx[s_blocks[0]];
    int end_block_idx   = block_id_to_idx[t_blocks[0]];

    std::vector<int> path_aps;
    std::vector<int> block_path_indices = findPathInBlockGraph(
        start_block_idx, end_block_idx, num_blocks, block_edges, path_aps);

    if (block_path_indices.empty()) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, 0, time);
    }

    // Build ordered block chain: block IDs and articulation points (entry/exit)
    // final_aps[0] = s, final_aps[k] = AP between block k-1 and k, final_aps[N] = t
    int chain_len = static_cast<int>(block_path_indices.size());
    std::vector<int> ordered_block_ids;
    for (int idx : block_path_indices)
        ordered_block_ids.push_back(block_nodes[idx].original_block_id);

    std::vector<int> final_aps(chain_len + 1, 0);
    final_aps[0] = s;
    for (int i = 0; i < static_cast<int>(path_aps.size()); ++i)
        final_aps[i + 1] = path_aps[i];
    final_aps[chain_len] = t;

    // Compute minimum diameter and reliability for each block in the chain
    std::vector<int> BlockDiam(chain_len, 0);
    int diamsum = 0;
    for (int i = 0; i < chain_len; ++i) {
        ReliabilityGraph block_graph;
        std::vector<int> map_new_to_orig, map_orig_to_new;
        getBlockGraphAndMap(ordered_block_ids[i], decomposition, block_graph,
                            map_new_to_orig, map_orig_to_new);

        int orig_entry = final_aps[i];
        int orig_exit  = final_aps[i + 1];
        if (orig_entry >= static_cast<int>(map_orig_to_new.size()) ||
            orig_exit  >= static_cast<int>(map_orig_to_new.size()) ||
            map_orig_to_new[orig_entry] == -1 || map_orig_to_new[orig_exit] == -1) {
            double time = (double)(clock() - start) / CLOCKS_PER_SEC;
            return ReliabilityResult(0.0, recursion_counter, time);
        }
        int ls = map_orig_to_new[orig_entry];
        int lt = map_orig_to_new[orig_exit];
        BlockDiam[i] = block_graph.calculateDistance(ls, lt, UpperBound);
        diamsum += BlockDiam[i];
    }

    if (diamsum > UpperBound) {
        double time = (double)(clock() - start) / CLOCKS_PER_SEC;
        return ReliabilityResult(0.0, recursion_counter, time);
    }

    int gap = UpperBound - diamsum;

    // Compute block reliabilities R_i(BlockDiam[i] + j) for j in [0, gap]
    std::vector<std::vector<double>> BlockReliab(chain_len, std::vector<double>(gap + 1, 0.0));
    for (int i = 0; i < chain_len; ++i) {
        ReliabilityGraph block_graph;
        std::vector<int> map_new_to_orig, map_orig_to_new;
        getBlockGraphAndMap(ordered_block_ids[i], decomposition, block_graph,
                            map_new_to_orig, map_orig_to_new);

        int ls = map_orig_to_new[final_aps[i]];
        int lt = map_orig_to_new[final_aps[i + 1]];

        for (int j = 0; j <= gap; ++j) {
            auto res = block_graph.calculateReliabilityBetweenVertices(ls, lt, BlockDiam[i] + j);
            BlockReliab[i][j] = res.reliability;
            recursion_counter += res.recursions;
        }
    }

    // Convert CDF → PMF for all blocks except the last
    for (int i = 0; i < chain_len - 1; ++i)
        for (int j = gap; j > 0; --j)
            BlockReliab[i][j] -= BlockReliab[i][j - 1];

    // Convolve block reliabilities (Migov's formula)
    std::vector<double> curBlockRel(gap + 1, 0.0);
    for (int i = chain_len - 1; i > 0; --i) {
        for (int diam = 0; diam <= gap; ++diam)
            for (int j = 0; j <= diam; ++j)
                curBlockRel[diam] += BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
        for (int j = 0; j <= gap; ++j) {
            BlockReliab[i - 1][j] = curBlockRel[j];
            curBlockRel[j] = 0.0;
        }
    }

    double reliability = BlockReliab[0][gap];
    double time = (double)(clock() - start) / CLOCKS_PER_SEC;

    LOG_INFO("Simple Factoring with Decomposition completed: reliability={}, recursions={}, time={}s",
             reliability, recursion_counter, time);
    return ReliabilityResult(reliability, recursion_counter, time);
}

ReliabilityResult ReliabilityGraph::calculateReliabilityWithParallelMDecomposition(VertexId s, VertexId t, int d) const {
    // Parallel version stub - alias to sequential for now
    LOG_DEBUG("Using parallel M-Decomposition (fallback to sequential)");
    return calculateReliabilityWithMDecomposition(s, t, d);
}

// -----------------------------------------------------------------------------
// Other Stubs/Helpers
// -----------------------------------------------------------------------------

ReliabilityGraph ReliabilityGraph::restoreKBlock(int b, const std::vector<int>& d) const {
    ReliabilityGraph g;
    std::vector<int> map1, map2;
    getBlockGraphAndMap(b, d, g, map1, map2);
    return g;
}

ReliabilityGraph ReliabilityGraph::restoreBlock(int b, const std::vector<int>& d) const {
    return restoreKBlock(b, d);
}

ReliabilityGraph ReliabilityGraph::createInducedSubgraph(const std::vector<VertexId>& v, int b) const {
    // Not needed for M-Decomp if we use getBlockGraphAndMap?
    return ReliabilityGraph();
}

bool ReliabilityGraph::areVerticesInSameBlock(const std::vector<VertexId>& v, int b, VertexId v1, VertexId v2) const {
    return false;
}

void ReliabilityGraph::outputToFile(std::ofstream& out) const {
    // CSV style output of graph structure?
    // Legacy kGraphFileOutput
}

void ReliabilityGraph::performFactoring(ReliabilityGraph g, int v, int d, double r) const {}
void ReliabilityGraph::perform2VertexFactoring(ReliabilityGraph g, VertexId s, VertexId t, int v, int d, double r) const {}
void ReliabilityGraph::performMFactoring(ReliabilityGraph g, VertexId s, VertexId t, int v, int d, double r, int l, int u) const {}
void ReliabilityGraph::performParallelMFactoring(ReliabilityGraph g, VertexId s, VertexId t, int v, int d, double r, int l, int u) const {}

} // namespace graph_reliability
