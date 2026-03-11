/**
 * @file Graph.h
 * @brief Base graph structure for reliability calculations
 * @author Graduate Work Project
 * @date 2026
 * 
 * This class provides the fundamental graph representation using Compressed Sparse Row (CSR)
 * format for efficient memory usage and fast adjacency queries. The graph stores:
 * - KAO: Array of offsets for each vertex in the adjacency list
 * - FO: Adjacency list (flat array of neighbors)
 * - PArray: Probability array for each edge
 * 
 * @note All vertex indices are 0-based internally
 * @note The graph is undirected (each edge appears twice in FO)
 */

#pragma once

#include <vector>
#include <memory>
#include <string>

namespace graph_reliability {

/**
 * @brief Base graph structure using Compressed Sparse Row (CSR) format
 * 
 * This class represents a graph using the CSR format for efficient memory usage
 * and fast adjacency queries. The graph stores:
 * - KAO: Array of offsets for each vertex in the adjacency list
 * - FO: Adjacency list (flat array of neighbors)
 * - PArray: Probability array for each edge
 */
class Graph {
public:
    using VertexId = int;
    using EdgeId = size_t;
    using Probability = double;

    // CSR format arrays
    std::vector<VertexId> kao_;      ///< Offset array for CSR format
    std::vector<VertexId> fo_;       ///< Adjacency list (flat array)
    std::vector<Probability> p_array_; ///< Edge probabilities

    /**
     * @brief Default constructor
     */
    Graph() = default;

    /**
     * @brief Constructor with CSR data
     * @param kao Offset array for CSR format
     * @param fo Adjacency list
     * @param p_array Edge probabilities
     */
    Graph(std::vector<VertexId> kao, 
          std::vector<VertexId> fo, 
          std::vector<Probability> p_array);

    /**
     * @brief Virtual destructor for inheritance
     */
    virtual ~Graph() = default;

    // Enable copy constructor and assignment operator
    Graph(const Graph&) = default;
    Graph& operator=(const Graph&) = default;

    // Enable move constructor and assignment operator
    Graph(Graph&&) = default;
    Graph& operator=(Graph&&) = default;

    /**
     * @brief Get number of vertices in the graph
     * @return Number of vertices
     */
    size_t numVertices() const noexcept;

    /**
     * @brief Get number of edges in the graph
     * @return Number of edges
     */
    size_t numEdges() const noexcept;

    /**
     * @brief Check if an edge exists between two vertices
     * @param from Source vertex
     * @param to Destination vertex
     * @return True if edge exists, false otherwise
     */
    bool hasEdge(VertexId from, VertexId to) const;

    /**
     * @brief Find edge index in the adjacency list
     * @param from Source vertex
     * @param to Destination vertex
     * @return Edge index or -1 if not found
     */
    EdgeId findEdge(VertexId from, VertexId to) const;

    /**
     * @brief Calculate shortest path distance using Dijkstra's algorithm
     * @param source Source vertex (0-based)
     * @param destination Destination vertex (0-based)
     * @param max_distance Maximum distance constant (early termination threshold)
     * @return Distance or max_distance+1 if no path exists or distance exceeds max_distance
     * 
     * @complexity O(E log V) where E is number of edges, V is number of vertices
     * @note Uses priority queue for efficient shortest path calculation
     * @note Early termination when distance exceeds max_distance for performance
     */
    int calculateDistance(VertexId source, VertexId destination, int max_distance) const;

    /**
     * @brief Remove an edge from the graph
     * @param from Source vertex
     * @param to Destination vertex
     * @return New graph without the specified edge
     */
    Graph removeEdge(VertexId from, VertexId to) const;

    /**
     * @brief Swap two vertices in the graph
     * @param vertex1 First vertex
     * @param vertex2 Second vertex
     * @return New graph with swapped vertices
     */
    Graph swapVertices(VertexId vertex1, VertexId vertex2) const;

    /**
     * @brief Decompose graph by removing specified vertices
     * @param vertices_to_remove Vertices to remove
     * @return Decomposition result vector
     */
    std::vector<int> decomposeByVertices(const std::vector<VertexId>& vertices_to_remove) const;

    /**
     * @brief Find articulation points in the graph
     * @param start_vertex Starting vertex for DFS
     * @param parent_vertex Parent vertex (-1 for root)
     * @return Vector of articulation points
     */
    std::vector<VertexId> findArticulationPoints(VertexId start_vertex, VertexId parent_vertex = -1) const;

    /**
     * @brief Search for minimum cut in the graph
     * @return Cut decomposition result
     */
    std::vector<int> findMinimumCut() const;

    /**
     * @brief Decompose graph into two components
     * @return Two-component decomposition result
     */
    std::vector<int> decomposeIntoTwoComponents() const;

    /**
     * @brief Convert edge list to KAO/FO format
     * @param input_path Path to input edge list file
     * @param output_path Path to output KAO/FO file
     * @param reliability Default reliability for all edges
     */
    virtual void convertEdgeListToKAOFO(const std::string& input_path,
                                       const std::string& output_path,
                                       Probability reliability = 0.9) const;

    /**
     * @brief Convert KAO/FO format to edge list
     * @param input_path Path to input KAO/FO file
     * @param output_path Path to output edge list file
     */
    virtual void convertKAOFOToEdgeList(const std::string& input_path,
                                       const std::string& output_path) const;

private:
    /**
     * @brief Validate graph structure
     * @throws std::invalid_argument if graph is invalid
     */
    void validateGraph() const;

    /**
     * @brief Normalize path separators
     * @param path Input path
     * @return Normalized path
     */
    static std::string normalizePath(const std::string& path);
};

} // namespace graph_reliability
