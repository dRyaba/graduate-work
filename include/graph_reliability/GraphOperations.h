/**
 * @file GraphOperations.h
 * @brief Graph operations and utilities
 * @author Graduate Work Project
 * @date 2026
 */

#pragma once

#include "Graph.h"
#include "ReliabilityGraph.h"
#include "DataImporter.h"
#include <vector>
#include <memory>

namespace graph_reliability {

/**
 * @brief Graph operations and utilities class
 * 
 * This class provides various graph operations and utility functions
 * for graph manipulation and analysis.
 */
class GraphOperations {
public:
    /**
     * @brief Merge two graphs at specified vertices
     * @param graph1 First graph
     * @param graph2 Second graph
     * @param connection_vertices Number of vertices to connect
     * @return Merged graph
     */
    static ReliabilityGraph mergeGraphs(const ReliabilityGraph& graph1,
                                       const ReliabilityGraph& graph2,
                                       int connection_vertices);

    /**
     * @brief Merge multiple graphs from file
     * @param connection_vertices Number of vertices to connect
     * @param data_importer Data importer instance
     * @return Merged graph
     */
    static ReliabilityGraph mergeGraphsFromFile(int connection_vertices,
                                               DataImporter& data_importer);

    /**
     * @brief Convert edge list to KAO/FO format
     * @param input_path Input edge list file path
     * @param output_path Output KAO/FO file path
     * @param reliability Default reliability for all edges
     */
    static void convertEdgeListToKAOFO(const std::string& input_path,
                                      const std::string& output_path,
                                      double reliability = 0.9);

    /**
     * @brief Convert KAO/FO format to edge list
     * @param input_path Input KAO/FO file path
     * @param output_path Output edge list file path
     */
    static void convertKAOFOToEdgeList(const std::string& input_path,
                                      const std::string& output_path);

private:
    /**
     * @brief Calculate new vertex mapping for graph merging
     * @param graph1_size Size of first graph
     * @param graph2_size Size of second graph
     * @param connection_vertices Number of connection vertices
     * @return Vertex mapping array
     */
    static std::vector<int> calculateVertexMapping(int graph1_size,
                                                  int graph2_size,
                                                  int connection_vertices);
};

} // namespace graph_reliability
