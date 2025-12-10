/**
 * @file TestSuite.cpp
 * @brief Implementation of the TestSuite class
 * @author Graduate Work Project
 * @date 2024
 */

#include "graph_reliability/TestSuite.h"
#include "graph_reliability/GraphOperations.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace graph_reliability {

TestSuite::TestSuite(DataImporter& data_importer)
    : data_importer_(data_importer), number_of_repetitions_(1) {
}

void TestSuite::setTestConfigurations(const std::map<std::string, std::vector<int>>& configs) {
    test_configurations_ = configs;
}

std::vector<std::string> TestSuite::getAvailableMethods() const {
    return {
        "Standard Factoring (Level 0)",
        "Recursive Decomposition (Level 1)",
        "Simple Factoring (Level 2)",
        "M-Decomposition (Level 3)"
    };
}

TestRunResult TestSuite::runSingleTest(const TestConfiguration& config, int method_id) {
    std::vector<double> times;
    std::vector<double> reliabilities;
    std::vector<long long> recursions;

    std::cout << "Testing " << config.graph_filename 
              << " (D=" << config.upper_bound_diameter << ")"
              << " Reps=" << config.number_of_repetitions << "..." << std::endl;

    for (int i = 0; i < config.number_of_repetitions; ++i) {
        try {
            // Load fresh graph for each run to ensure clean state
            auto graph = data_importer_.loadKAOGraph(config.graph_filename);
            
            if (!graph) {
                std::cerr << "Failed to load graph: " << config.graph_filename << std::endl;
                continue;
            }

            int s = config.source_vertex;
            int t = config.target_vertex;

            ReliabilityResult result;
            
            switch (method_id) {
                case 0: // Level 0: Standard Factoring (Baseline - SLOWEST)
                    result = graph->calculateReliabilityBetweenVertices(s, t, config.upper_bound_diameter);
                    break;
                case 1: // Level 1: Recursive Decomposition (Nested Recursion - INEFFICIENT)
                    result = graph->calculateReliabilityWithRecursiveDecomposition(s, t, config.upper_bound_diameter);
                    break;
                case 2: // Level 2: Simple Factoring (Convolution - FAST)
                    result = graph->calculateReliabilityWithDecomposition(s, t, config.upper_bound_diameter);
                    break;
                case 3: // Level 3: M-Decomposition (FASTEST)
                    result = graph->calculateReliabilityWithMDecomposition(s, t, config.upper_bound_diameter);
                    break;
                default:
                    std::cerr << "Unknown method ID: " << method_id << std::endl;
                    return TestRunResult();
            }

            times.push_back(result.execution_time_sec);
            reliabilities.push_back(result.reliability);
            recursions.push_back(result.recursions);

        } catch (const std::exception& e) {
            std::cerr << "Error in test run " << i + 1 << ": " << e.what() << std::endl;
        }
    }

    if (times.empty()) return TestRunResult();

    double avg_time = 0;
    double avg_rel = 0;
    long long avg_rec = 0;

    for (double t : times) avg_time += t;
    for (double r : reliabilities) avg_rel += r;
    for (long long rec : recursions) avg_rec += rec;

    return TestRunResult(
        avg_time / times.size(),
        avg_rel / reliabilities.size(),
        avg_rec / recursions.size()
    );
}

bool TestSuite::runComprehensiveTests(int method_id, const std::string& output_filename) {
    std::vector<std::pair<TestConfiguration, TestRunResult>> results;
    bool all_success = true;
    std::vector<std::string> methods = getAvailableMethods();
    std::string method_name = (method_id >= 0 && method_id < methods.size()) ? methods[method_id] : "Unknown";

    for (const auto& [filename, diameters] : test_configurations_) {
        try {
            // Find s/t for this graph
            auto temp_graph = data_importer_.loadKAOGraph(filename);
            if (!temp_graph) {
                std::cerr << "Skipping " << filename << ": could not load" << std::endl;
                all_success = false;
                continue;
            }

            auto [s, t] = findSourceAndTargetVertices(*temp_graph);
            if (s == -1 || t == -1) {
                std::cerr << "Skipping " << filename << ": could not find source/target vertices" << std::endl;
                all_success = false;
                continue;
            }

            for (int d : diameters) {
                TestConfiguration config(filename, s, t, d, number_of_repetitions_);
                TestRunResult result = runSingleTest(config, method_id);
                results.push_back({config, result});
            }

        } catch (const std::exception& e) {
            std::cerr << "Error processing " << filename << ": " << e.what() << std::endl;
            all_success = false;
        }
    }

    writeResultsToCSV(results, method_name, output_filename);
    return all_success;
}

std::pair<int, int> TestSuite::findSourceAndTargetVertices(const ReliabilityGraph& graph) const {
    // Logic from legacy main.cpp:
    // Find vertices where Targets[i] == 1. First is s, second is t.
    // Note: ReliabilityGraph uses 0-based target_vertices_ vector.
    // But legacy Targets was 1-based vector (size N+1).
    // My DataImporter converts Targets to 0-based vector of size N?
    // Or size N+1? 
    // Let's check DataImporter.cpp: 
    // "if (!targets.empty()) targets.erase(targets.begin());" -> Size N.
    // So indices are 0..N-1.

    int s = -1, t = -1;
    const auto& targets = graph.target_vertices_;
    
    for (size_t i = 0; i < targets.size(); ++i) {
        if (targets[i] == 1) {
            if (s == -1) s = static_cast<int>(i);
            else if (t == -1) t = static_cast<int>(i);
        }
    }
    
    return {s, t};
}

void TestSuite::writeResultsToCSV(const std::vector<std::pair<TestConfiguration, TestRunResult>>& results,
                                 const std::string& method_name,
                                 const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Could not open output file: " << filename << std::endl;
        return;
    }

    file << "Graph,S,T,D,Reps,Method,AvgTime,AvgRel,AvgRecs\n";
    
    for (const auto& [config, result] : results) {
        file << config.graph_filename << ","
             << config.source_vertex << ","
             << config.target_vertex << ","
             << config.upper_bound_diameter << ","
             << config.number_of_repetitions << ","
             << method_name << ","
             << std::fixed << std::setprecision(6) << result.average_time_seconds << ","
             << std::setprecision(9) << result.average_reliability << ","
             << result.average_recursions << "\n";
    }
    
    std::cout << "Results written to " << filename << std::endl;
}

} // namespace graph_reliability

