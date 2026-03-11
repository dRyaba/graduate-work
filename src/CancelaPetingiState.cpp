/**
 * @file CancelaPetingiState.cpp
 * @brief Implementation of CPFM state
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/CancelaPetingiState.h"
#include "graph_reliability/Graph.h"
#include <algorithm>
#include <set>

namespace graph_reliability {

size_t CancelaPetingiState::countFeasiblePathsInP(EdgeId e) const {
    auto it = P_of_edge_.find(e);
    if (it == P_of_edge_.end()) return 0;
    size_t count = 0;
    for (size_t pi : it->second) {
        if (pi < feasible_.size() && feasible_[pi]) count++;
    }
    return count;
}

void CancelaPetingiState::buildPOfEdge(const std::vector<Path>& paths) {
    P_of_edge_.clear();
    for (size_t pi = 0; pi < paths.size(); ++pi) {
        for (EdgeId e : paths[pi]) {
            P_of_edge_[e].push_back(pi);
        }
    }
}

CancelaPetingiState::CancelaPetingiState(const std::vector<Path>& paths, const Graph& graph)
    : paths_(paths) {
    linksp_.resize(paths.size());
    feasible_.resize(paths.size(), true);
    npst_ = static_cast<int>(paths.size());
    connected_st_ = false;

    std::set<EdgeId> edge_set;
    for (const auto& path : paths) {
        for (EdgeId e : path) {
            edge_set.insert(e);
        }
    }

    for (EdgeId e : edge_set) {
        if (e >= graph.p_array_.size()) continue;
        double p = graph.p_array_[e];
        edge_reliability_[e] = p;
        if (p > 0 && p < 1.0) {
            active_edges_.push_back(e);
        }
    }

    for (size_t pi = 0; pi < paths.size(); ++pi) {
        int count = 0;
        for (EdgeId e : paths[pi]) {
            auto it = edge_reliability_.find(e);
            if (it != edge_reliability_.end() && it->second < 1.0) {
                count++;
            }
        }
        linksp_[pi] = count;
    }

    buildPOfEdge(paths);
}

void CancelaPetingiState::applyContract(EdgeId e) {
    auto it = P_of_edge_.find(e);
    if (it == P_of_edge_.end()) return;

    for (size_t pi : it->second) {
        if (!feasible_[pi]) continue;
        linksp_[pi]--;
        if (!connected_st_ && linksp_[pi] == 0) {
            connected_st_ = true;
        }
    }
}

void CancelaPetingiState::applyDelete(EdgeId e) {
    auto it = P_of_edge_.find(e);
    if (it == P_of_edge_.end()) return;

    for (size_t pi : it->second) {
        if (feasible_[pi]) {
            feasible_[pi] = false;
            npst_--;
        }
    }
}

std::vector<CancelaPetingiState::EdgeId> CancelaPetingiState::getActiveEdgesByPCount() const {
    std::vector<EdgeId> result;
    for (const auto& [e, r] : edge_reliability_) {
        if (r > 0 && r < 1.0) {
            result.push_back(e);
        }
    }
    std::sort(result.begin(), result.end(), [this](EdgeId a, EdgeId b) {
        size_t pa = countFeasiblePathsInP(a);
        size_t pb = countFeasiblePathsInP(b);
        if (pa != pb) return pa > pb;
        return a < b;
    });
    return result;
}

CancelaPetingiState::EdgeId CancelaPetingiState::selectPivotEdge() const {
    auto edges = getActiveEdgesByPCount();
    return edges.empty() ? static_cast<EdgeId>(-1) : edges.front();
}

int CancelaPetingiState::applyISPT(EdgeId e) {
    auto it_e = P_of_edge_.find(e);
    if (it_e == P_of_edge_.end()) return 0;

    const std::vector<size_t>& Pe = it_e->second;
    
    std::set<size_t> feasible_Pe;
    for (size_t pi : Pe) {
        if (feasible_[pi]) feasible_Pe.insert(pi);
    }
    if (feasible_Pe.empty()) return 0;
    
    int merged = 0;

    std::vector<EdgeId> to_merge;
    for (const auto& [f, Pf] : P_of_edge_) {
        if (f == e) continue;
        auto rit = edge_reliability_.find(f);
        if (rit == edge_reliability_.end() || rit->second <= 0 || rit->second >= 1.0) continue;

        std::set<size_t> feasible_Pf;
        for (size_t pi : Pf) {
            if (feasible_[pi]) feasible_Pf.insert(pi);
        }
        
        if (feasible_Pf != feasible_Pe) continue;
        
        to_merge.push_back(f);
    }

    for (EdgeId f : to_merge) {
        auto rit = edge_reliability_.find(f);
        if (rit == edge_reliability_.end()) continue;
        double rf = rit->second;
        edge_reliability_[e] *= rf;
        edge_reliability_.erase(f);
        P_of_edge_.erase(f);
        for (size_t pi : feasible_Pe) {
            linksp_[pi]--;
        }
        merged++;
    }

    return merged;
}

int CancelaPetingiState::applyGlobalISPT() {
    int total_merged = 0;
    
    std::vector<EdgeId> active;
    for (const auto& [e, r] : edge_reliability_) {
        if (r > 0 && r < 1.0) {
            active.push_back(e);
        }
    }
    
    std::map<EdgeId, std::set<size_t>> feasible_P;
    for (EdgeId e : active) {
        auto it = P_of_edge_.find(e);
        if (it == P_of_edge_.end()) continue;
        for (size_t pi : it->second) {
            if (feasible_[pi]) feasible_P[e].insert(pi);
        }
    }
    
    std::set<EdgeId> already_merged;
    
    for (size_t i = 0; i < active.size(); ++i) {
        EdgeId e = active[i];
        if (already_merged.count(e)) continue;
        if (feasible_P[e].empty()) continue;
        
        std::vector<EdgeId> to_merge;
        for (size_t j = i + 1; j < active.size(); ++j) {
            EdgeId f = active[j];
            if (already_merged.count(f)) continue;
            
            if (feasible_P[f] == feasible_P[e]) {
                to_merge.push_back(f);
            }
        }
        
        for (EdgeId f : to_merge) {
            auto rit = edge_reliability_.find(f);
            if (rit == edge_reliability_.end()) continue;
            
            double rf = rit->second;
            edge_reliability_[e] *= rf;
            edge_reliability_.erase(f);
            P_of_edge_.erase(f);
            
            for (size_t pi : feasible_P[e]) {
                linksp_[pi]--;
            }
            
            already_merged.insert(f);
            total_merged++;
        }
    }
    
    return total_merged;
}

std::vector<CancelaPetingiState::EdgeId> CancelaPetingiState::getSortedEdgesForFactoring() const {
    std::vector<std::pair<EdgeId, size_t>> edges_with_count;
    
    for (const auto& [e, r] : edge_reliability_) {
        if (r > 0 && r < 1.0) {
            auto it = P_of_edge_.find(e);
            size_t count = (it != P_of_edge_.end()) ? it->second.size() : 0;
            edges_with_count.emplace_back(e, count);
        }
    }
    
    std::sort(edges_with_count.begin(), edges_with_count.end(),
              [](const auto& a, const auto& b) {
                  if (a.second != b.second) return a.second > b.second;
                  return a.first < b.first;
              });
    
    std::vector<EdgeId> result;
    result.reserve(edges_with_count.size());
    for (const auto& [e, _] : edges_with_count) {
        result.push_back(e);
    }
    
    return result;
}

} // namespace graph_reliability
