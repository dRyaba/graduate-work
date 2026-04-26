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
#include <functional>
#include <optional>

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
 * @brief Status of a single cross-check cell.
 */
enum class CrossCheckStatus { OK, TIMEOUT, ERROR };

/**
 * @brief One row in the cross-check CSV: (graph, s, t, d, method) → result or timeout/error.
 */
struct CrossCheckRow {
    std::string graph;
    int s = 0;
    int t = 0;
    int diameter = 0;
    int method_id = 0;
    std::string method_name;
    CrossCheckStatus status = CrossCheckStatus::OK;
    double reliability = 0.0;
    double time_seconds = 0.0;
    long long recursions = 0;
    std::string error_message;
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
     * @brief Run cross-consistency check across all methods m0..m5.
     *
     * For each (graph, diameter) in the current test_configurations_ (plus
     * optional real-world networks), runs every method 0..5 with a per-cell
     * timeout and compares the reliabilities against a baseline (first
     * successful method in order m0, m3, m5). Writes a detailed CSV and
     * prints a summary. Returns true if no discrepancies exceed `tolerance`
     * (timeouts are not counted as failures).
     *
     * @param timeout_sec Per-cell timeout in seconds
     * @param output_filename Output CSV path
     * @param include_real_networks If true, also test Geant2004 / IEEE-118
     * @param tolerance Max allowed |R - R_baseline| to consider methods consistent
     */
    bool runCrossCheck(int timeout_sec,
                       const std::string& output_filename,
                       bool include_real_networks,
                       double tolerance = 1e-10,
                       int d_count = 4,
                       int d_step = 1,
                       const std::vector<int>& active_methods = {0, 1, 2, 3, 4, 5});

    /**
     * @brief Choose content-rich diameters for cell (s, t).
     *
     * Returns { dist(s,t), dist+step, dist+2*step, ..., dist+(count-1)*step },
     * clipped to [1, |V|-1]. If s and t are unreachable within |V|-1 hops,
     * returns an empty vector. Degenerate d < dist(s,t) is never tested,
     * so every cell is guaranteed R > 0 (until saturation).
     *
     * @param graph Graph to measure dist on
     * @param s Source vertex (0-based)
     * @param t Target vertex (0-based)
     * @param count Number of diameters to produce (>= 1)
     * @param step Spacing between successive diameters (>= 1)
     * @return Sorted strictly-increasing diameters, or empty on unreachable.
     */
    static std::vector<int> chooseDiameters(const ReliabilityGraph& graph,
                                            int s, int t,
                                            int count = 4, int step = 1);

    /**
     * @brief Run a reliability calculation in a detached worker thread with a timeout.
     *
     * If the worker finishes in time, returns its ReliabilityResult. If it does
     * not, returns std::nullopt; the worker thread is detached and will be torn
     * down when the process exits. This is acceptable for one-shot cross-check
     * runs and CI tests where adding cooperative cancellation to every algorithm
     * would be disproportionate. Exceptions inside the worker are swallowed and
     * reported as errors via out_error.
     *
     * @param calc Nullary callable returning ReliabilityResult
     * @param timeout_sec Maximum wait in seconds (>= 1)
     * @param out_error Optional out-param: set to exception message on ERROR
     * @return ReliabilityResult on success, nullopt on timeout or exception
     */
    static std::optional<ReliabilityResult> runWithTimeout(
        std::function<ReliabilityResult()> calc,
        int timeout_sec,
        std::string* out_error = nullptr);

    /**
     * @brief Multi-diameter twin of runWithTimeout for CDF-returning methods.
     *
     * Same contract: runs `calc` on a worker thread, returns its
     * `ReliabilityCdfResult` if it finishes inside `timeout_sec`, otherwise
     * detaches the worker and returns nullopt.
     */
    static std::optional<ReliabilityCdfResult> runWithTimeoutCdf(
        std::function<ReliabilityCdfResult()> calc,
        int timeout_sec,
        std::string* out_error = nullptr);

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

    /**
     * @brief Write a vector of cross-check rows to CSV.
     */
    static void writeCrossCheckCSV(const std::vector<CrossCheckRow>& rows,
                                   const std::string& filename);
};

} // namespace graph_reliability
