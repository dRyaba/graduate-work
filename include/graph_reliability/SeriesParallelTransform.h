/**
 * @file SeriesParallelTransform.h
 * @brief Series-parallel transformation for diameter-constrained reliability
 * @author Graduate Work Project
 * @date 2026
 *
 * Removes chains (2-degree non-terminal nodes) before path enumeration.
 * For DCR, parallel transformation is not applicable (NesterovMigov).
 *
 * @see NesterovMigov.pdf Section III-B
 */

#pragma once

#include "Graph.h"
#include <utility>

namespace graph_reliability {

/**
 * @brief Apply series transformation to reduce chains
 *
 * Removes 2-degree non-terminal vertices, replacing (s,v)-(v,t) with (s,t)
 * when r(s,t) = r(s,v) * r(v,t). Does not create multi-edges (parallel not applicable for DCR).
 *
 * @param graph Input graph
 * @param terminals Vertices that are terminals (not reduced)
 * @return Reduced graph and multiplicative factor (product of removed edge probs for chains we couldn't merge)
 */
std::pair<Graph, double> applySeriesTransformation(const Graph& graph,
                                                   const std::vector<int>& terminals);

} // namespace graph_reliability
