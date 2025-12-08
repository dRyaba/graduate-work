/**
 * @file GraphOperations.cpp
 * @brief Implementation of GraphOperations class
 * @author Graduate Work Project
 * @date 2024
 */

#include "graph_reliability/GraphOperations.h"
#include <iostream>

namespace graph_reliability {

ReliabilityGraph GraphOperations::mergeGraphs(const ReliabilityGraph& graph1,
                                             const ReliabilityGraph& graph2,
                                             int connection_vertices) {
    // Stub for now
    throw std::runtime_error("mergeGraphs not implemented yet");
}

ReliabilityGraph GraphOperations::mergeGraphsFromFile(int connection_vertices,
                                                     DataImporter& data_importer) {
    auto graphs = data_importer.loadGraphsToMerge();
    if (graphs.size() < 2) {
        throw std::runtime_error("Need at least 2 graphs to merge");
    }
    // return mergeGraphs(*graphs[0], *graphs[1], connection_vertices);
    throw std::runtime_error("mergeGraphsFromFile not implemented yet");
}

void GraphOperations::convertEdgeListToKAOFO(const std::string& input_path,
                                            const std::string& output_path,
                                            double reliability) {
    Graph g;
    g.convertEdgeListToKAOFO(input_path, output_path, reliability);
}

void GraphOperations::convertKAOFOToEdgeList(const std::string& input_path,
                                            const std::string& output_path) {
    Graph g;
    g.convertKAOFOToEdgeList(input_path, output_path);
}

std::vector<int> GraphOperations::calculateVertexMapping(int graph1_size,
                                                        int graph2_size,
                                                        int connection_vertices) {
    return {};
}

} // namespace graph_reliability

