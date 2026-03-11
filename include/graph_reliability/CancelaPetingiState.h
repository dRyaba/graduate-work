/**
 * @file CancelaPetingiState.h
 * @brief State structure for Cancela-Petingi path-based factoring
 * @author Graduate Work Project
 * @date 2026
 *
 * Holds path lists, P(e) index, and feasibility flags for CPFM recursion.
 *
 * @see TR0103.pdf (Cancela, Petingi), NesterovMigov.pdf
 */

#pragma once

#include "PathEnumerator.h"
#include <map>
#include <vector>
#include <cstddef>

namespace graph_reliability {

/**
 * @brief State for Cancela-Petingi factoring method (2-terminal)
 *
 * Tracks paths, their feasibility, and connectivity for recursive factoring.
 */
class CancelaPetingiState {
public:
    using EdgeId = PathEnumerator::EdgeId;
    using Path = PathEnumerator::Path;

    /** Paths Pst(d) - list of paths between s and t */
    std::vector<Path> paths_;

    /** P(e): for each edge e, indices of paths in paths_ that contain e */
    std::map<EdgeId, std::vector<size_t>> P_of_edge_;

    /** linksp: for each path, number of non-perfect (unreliable) edges */
    std::vector<int> linksp_;

    /** feasiblep: for each path, true if path has no failed edges */
    std::vector<bool> feasible_;

    /** npst: number of feasible paths between s and t */
    int npst_ = 0;

    /** connectedst: s and t connected by a perfect path */
    bool connected_st_ = false;

    /** edge_reliability_: current reliability for each edge (changes with SPT merge) */
    std::map<EdgeId, double> edge_reliability_;

    /** Edges that are still active (0 < re < 1) */
    std::vector<EdgeId> active_edges_;

    /**
     * @brief Build state from enumerated paths and graph edge probabilities
     * @param paths Paths from PathEnumerator::enumeratePaths
     * @param graph Graph for edge probabilities and structure
     */
    CancelaPetingiState(const std::vector<Path>& paths, const Graph& graph);

    /** Default copy and move */
    CancelaPetingiState() = default;
    CancelaPetingiState(const CancelaPetingiState&) = default;
    CancelaPetingiState& operator=(const CancelaPetingiState&) = default;
    CancelaPetingiState(CancelaPetingiState&&) = default;
    CancelaPetingiState& operator=(CancelaPetingiState&&) = default;

    /**
     * @brief Apply Contract branch: edge e is reliable
     * Updates linksp, connected_st for paths in P(e)
     */
    void applyContract(EdgeId e);

    /**
     * @brief Apply Delete branch: edge e is failed
     * Updates feasible_, npst_ for paths in P(e)
     */
    void applyDelete(EdgeId e);

    /** @brief All terminals connected by perfect paths (2-terminal: connected_st_) */
    bool isTerminalSuccess() const { return connected_st_; }

    /** @brief No feasible path remains (npst == 0) */
    bool isTerminalFail() const { return npst_ <= 0; }

    /**
     * @brief Get edges with 0 < re < 1, sorted by |P(e)| descending (ESS)
     */
    std::vector<EdgeId> getActiveEdgesByPCount() const;

    /**
     * @brief Select pivot edge: max |P(e)| among active edges
     * @return EdgeId or invalid if no active edge
     */
    EdgeId selectPivotEdge() const;

    /**
     * @brief Apply ISPT: merge edges f where P(f) == P(e) into e
     * Updates edge_reliability_[e] *= rf, linksp for paths in P(e)
     * @return Number of edges merged
     */
    int applyISPT(EdgeId e);

    /**
     * @brief Apply global ISPT: merge ALL edges with identical P(e) sets
     * Should be called once before recursion (Nesterov's newReduce approach)
     * More aggressive than per-pivot ISPT
     * @return Total number of edges merged
     */
    int applyGlobalISPT();

    /**
     * @brief Get edges sorted by |P(e)| descending (fixed order for ESS)
     * Call once before recursion, then iterate by index
     * @return Sorted edge list (fixed order)
     */
    std::vector<EdgeId> getSortedEdgesForFactoring() const;

private:
    void buildPOfEdge(const std::vector<Path>& paths);

    /** @brief Count feasible paths in P(e) for ESS */
    size_t countFeasiblePathsInP(EdgeId e) const;
};

} // namespace graph_reliability
