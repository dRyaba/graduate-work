/**
 * @file TestSuite.h
 * @brief Test suite for graph reliability calculations
 * @author Graduate Work Project
 * @date 2026
 */

#pragma once

#include "ReliabilityGraph.h"
#include "DataImporter.h"
#include <string>
#include <vector>
#include <map>
#include <chrono>

namespace graph_reliability {

/**
 * @brief Test configuration structure
 */
struct TestConfiguration {
    std::string graph_filename;        ///< Graph file name
    int source_vertex;                 ///< Source vertex (1-based)
    int target_vertex;                 ///< Target vertex (1-based)
    int upper_bound_diameter;          ///< Upper bound for diameter
    int number_of_repetitions;         ///< Number of test repetitions

    TestConfiguration(const std::string& filename = "",
                     int source = 0,
                     int target = 0,
                     int diameter = 0,
                     int repetitions = 1)
        : graph_filename(filename), source_vertex(source), target_vertex(target),
          upper_bound_diameter(diameter), number_of_repetitions(repetitions) {}
};

/**
 * @brief Test run result structure
 */
struct TestRunResult {
    double average_time_seconds;       ///< Average execution time
    double average_reliability;        ///< Average reliability value
    long long average_recursions;      ///< Average number of recursions

    TestRunResult(double time = 0.0, double reliability = 0.0, long long recursions = 0)
        : average_time_seconds(time), average_reliability(reliability), 
          average_recursions(recursions) {}
};

/**
 * @brief Test suite for graph reliability calculations
 * 
 * This class provides a comprehensive testing framework for different
 * reliability calculation methods and algorithms.
 */
class TestSuite {
public:
    /**
     * @brief Constructor
     * @param data_importer Reference to data importer
     */
    explicit TestSuite(DataImporter& data_importer);
    
    /**
     * @brief Destructor
     */
    ~TestSuite() = default;
    
    // Disable copy constructor and assignment operator
    TestSuite(const TestSuite&) = delete;
    TestSuite& operator=(const TestSuite&) = delete;
    
    // Enable move constructor, disable move assignment (due to reference member)
    TestSuite(TestSuite&&) = default;
    TestSuite& operator=(TestSuite&&) = delete;

    /**
     * @brief Run single test configuration
     * @param config Test configuration
     * @param method_id Method ID (0=MDecompose, 1=Recursive, 2=Simple, 3=SF)
     * @return Test run result
     */
    TestRunResult runSingleTest(const TestConfiguration& config, int method_id);

    /**
     * @brief Run comprehensive test suite
     * @param method_id Method ID to test
     * @param output_filename Output CSV filename
     * @return True if all tests passed successfully
     */
    bool runComprehensiveTests(int method_id, const std::string& output_filename = "test_results.csv");

    /**
     * @brief Set test configurations
     * @param configs Map of filename to diameter bounds
     */
    void setTestConfigurations(const std::map<std::string, std::vector<int>>& configs);

    /**
     * @brief Get available test methods
     * @return Vector of method names
     */
    std::vector<std::string> getAvailableMethods() const;

    /**
     * @brief Set number of repetitions for tests
     * @param repetitions Number of repetitions
     */
    void setNumberOfRepetitions(int repetitions) noexcept { number_of_repetitions_ = repetitions; }

    /**
     * @brief Get number of repetitions
     * @return Number of repetitions
     */
    int getNumberOfRepetitions() const noexcept { return number_of_repetitions_; }

private:
    DataImporter& data_importer_;
    std::map<std::string, std::vector<int>> test_configurations_;
    int number_of_repetitions_;

    /**
     * @brief Find source and target vertices in graph
     * @param graph Graph to analyze
     * @return Pair of (source_vertex, target_vertex)
     */
    std::pair<int, int> findSourceAndTargetVertices(const ReliabilityGraph& graph) const;

    /**
     * @brief Execute reliability calculation method
     * @param graph Graph to test
     * @param source_vertex Source vertex
     * @param target_vertex Target vertex
     * @param diameter Upper bound diameter
     * @param method_id Method ID
     * @return Reliability result
     */
    ReliabilityResult executeMethod(const ReliabilityGraph& graph,
                                   int source_vertex,
                                   int target_vertex,
                                   int diameter,
                                   int method_id) const;

    /**
     * @brief Write test results to CSV file
     * @param results Test results
     * @param method_name Method name
     * @param filename Output filename
     */
    void writeResultsToCSV(const std::vector<std::pair<TestConfiguration, TestRunResult>>& results,
                          const std::string& method_name,
                          const std::string& filename) const;
};

} // namespace graph_reliability
