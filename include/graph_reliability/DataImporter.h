/**
 * @file DataImporter.h
 * @brief Centralized data import system for graph reliability calculations
 * @author Graduate Work Project
 * @date 2026
 */

#pragma once

#include "ReliabilityGraph.h"
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

namespace graph_reliability {

/**
 * @brief Centralized data import system for the graph reliability project
 * 
 * This class provides a unified interface for loading graphs from various formats
 * and managing data file paths. It replaces the old hardcoded path system with
 * a flexible, cross-platform solution.
 */
class DataImporter {
public:
    /**
     * @brief Constructor with base data path configuration
     * @param base_data_path Base path to data directory (default: "graphs_data/")
     */
    explicit DataImporter(const std::string& base_data_path = "graphs_data/");
    
    /**
     * @brief Destructor
     */
    ~DataImporter() = default;
    
    // Disable copy constructor and assignment operator
    DataImporter(const DataImporter&) = delete;
    DataImporter& operator=(const DataImporter&) = delete;
    
    // Enable move constructor and assignment operator
    DataImporter(DataImporter&&) = default;
    DataImporter& operator=(DataImporter&&) = default;
    
    /**
     * @brief Load graph from KAO format
     * @param filename File name (without path)
     * @return Smart pointer to loaded graph
     * @throws std::runtime_error if file not found or corrupted
     */
    std::unique_ptr<ReliabilityGraph> loadKAOGraph(const std::string& filename);
    
    /**
     * @brief Load graph from Edge List format
     * @param filename File name (without path)
     * @param reliability Reliability value for all edges (default: 0.9)
     * @return Smart pointer to loaded graph
     * @throws std::runtime_error if file not found or corrupted
     */
    std::unique_ptr<ReliabilityGraph> loadEdgeListGraph(const std::string& filename, 
                                                        double reliability = 0.9);
    
    /**
     * @brief Load multiple graphs from GraphsToMerge.txt file
     * @return Vector of smart pointers to loaded graphs
     * @throws std::runtime_error if file not found or corrupted
     */
    std::vector<std::unique_ptr<ReliabilityGraph>> loadGraphsToMerge();
    
    /**
     * @brief Get full path to data file
     * @param filename File name
     * @return Full path to file
     */
    std::string getFullPath(const std::string& filename) const;
    
    /**
     * @brief Check if data file exists
     * @param filename File name
     * @return True if file exists, false otherwise
     */
    bool fileExists(const std::string& filename) const;
    
    /**
     * @brief Get list of all available data files
     * @return Vector of file names
     */
    std::vector<std::string> getAvailableFiles() const;
    
    /**
     * @brief Set new base data path
     * @param new_base_path New base path
     */
    void setBasePath(const std::string& new_base_path);
    
    /**
     * @brief Get current base path
     * @return Current base path
     */
    std::string getBasePath() const noexcept { return base_data_path_; }
    
    /**
     * @brief Convert Edge List to KAO format and save
     * @param input_filename Input file name (Edge List)
     * @param output_filename Output file name (KAO)
     * @param reliability Reliability value for all edges
     * @throws std::runtime_error if conversion fails
     */
    void convertEdgeListToKAO(const std::string& input_filename, 
                             const std::string& output_filename, 
                             double reliability = 0.9);
    
    /**
     * @brief Convert KAO to Edge List format and save
     * @param input_filename Input file name (KAO)
     * @param output_filename Output file name (Edge List)
     * @throws std::runtime_error if conversion fails
     */
    void convertKAOToEdgeList(const std::string& input_filename, 
                             const std::string& output_filename);

private:
    std::string base_data_path_;
    
    /**
     * @brief Internal function to load KAO graph from stream
     * @param input_stream Input stream
     * @return Smart pointer to loaded graph
     * @throws std::runtime_error if data is corrupted
     */
    std::unique_ptr<ReliabilityGraph> loadKAOFromStream(std::ifstream& input_stream);
    
    /**
     * @brief Validate loaded graph
     * @param graph Pointer to graph for validation
     * @param filename File name (for error messages)
     * @throws std::runtime_error if graph is invalid
     */
    void validateGraph(const ReliabilityGraph* graph, const std::string& filename) const;
    
    /**
     * @brief Normalize path (replace separators, remove extra characters)
     * @param path Path to normalize
     * @return Normalized path
     */
    std::string normalizePath(const std::string& path) const;
};

} // namespace graph_reliability
