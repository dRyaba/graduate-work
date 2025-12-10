/**
 * @file main.cpp
 * @brief Main application for Graph Reliability Analysis
 * @author Graduate Work Project
 * @date 2024
 * 
 * This is the main entry point for the graph reliability analysis application.
 * It provides command-line interface for running reliability tests and
 * graph format conversions.
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include "graph_reliability.h"

using namespace graph_reliability;

/**
 * @brief Print application usage information
 */
void printUsage() {
    std::cout << "Graph Reliability Analysis Tool v" << getVersion() << "\n\n";
    std::cout << "Usage:\n";
    std::cout << "  ./graph_reliability [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --run <file> <s> <t> <d> <method> [reps]           Run test on specific graph file\n";
    std::cout << "  --test <method_id> [output_file]                   Run comprehensive tests on predefined graphs\n";
    std::cout << "  --convert edge2kao <input> <output> <reliability>  Convert Edge List to KAO format\n";
    std::cout << "  --convert kao2edge <input> <output>                Convert KAO format to Edge List\n";
    std::cout << "  --help                                             Show this help message\n\n";
    std::cout << "Test Methods (aligned with Optimization Levels):\n";
    std::cout << "  0 - Standard Factoring      (Level 0: Global, No Decomposition - SLOWEST)\n";
    std::cout << "  1 - Recursive Decomposition (Level 1: Nested Recursion - INEFFICIENT)\n";
    std::cout << "  2 - Simple Factoring        (Level 2: Convolution + Simple Facto - FAST)\n";
    std::cout << "  3 - M-Decomposition         (Level 3: Convolution + Modified Facto - FASTEST)\n\n";
    std::cout << "Parameters for --run:\n";
    std::cout << "  <file>    - Graph file path (KAO format)\n";
    std::cout << "  <s>       - Source vertex (0-based index, or -1 to auto-detect from targets)\n";
    std::cout << "  <t>       - Target vertex (0-based index, or -1 to auto-detect from targets)\n";
    std::cout << "  <d>       - Diameter upper bound\n";
    std::cout << "  <method>  - Method ID (0-3)\n";
    std::cout << "  [reps]    - Number of repetitions (default: 1)\n\n";
    std::cout << "Examples:\n";
    std::cout << "  ./graph_reliability --run graphs_data/K4_kao.txt 0 3 2 3\n";
    std::cout << "  ./graph_reliability --run graphs_data/3_blocks_sausage_3x3_kao.txt -1 -1 10 3 5\n";
    std::cout << "  ./graph_reliability --test 3 results.csv\n";
    std::cout << "  ./graph_reliability --convert edge2kao input.edgelist output.kao 0.9\n";
}

/**
 * @brief Handle single graph test run
 * @param args Command line arguments
 * @return Exit code
 */
int handleRun(const std::vector<std::string>& args) {
    // --run <file> <s> <t> <d> <method> [reps]
    if (args.size() < 6) {
        std::cerr << "Error: Insufficient arguments for --run\n";
        std::cerr << "Usage: --run <file> <s> <t> <d> <method> [reps]\n";
        return 1;
    }

    const std::string& filename = args[1];
    int s_vertex, t_vertex, diameter, method_id;
    int repetitions = 1;

    try {
        s_vertex = std::stoi(args[2]);
        t_vertex = std::stoi(args[3]);
        diameter = std::stoi(args[4]);
        method_id = std::stoi(args[5]);
        
        if (args.size() > 6) {
            repetitions = std::stoi(args[6]);
        }

        if (method_id < 0 || method_id > 3) {
            std::cerr << "Error: Method ID must be between 0 and 3\n";
            return 1;
        }
        if (diameter < 1) {
            std::cerr << "Error: Diameter must be positive\n";
            return 1;
        }
        if (repetitions < 1) {
            std::cerr << "Error: Repetitions must be positive\n";
            return 1;
        }
    } catch (const std::exception&) {
        std::cerr << "Error: Invalid numeric argument\n";
        return 1;
    }

    try {
        // Determine base path from filename
        std::string base_path = "";
        size_t last_sep = filename.find_last_of("/\\");
        if (last_sep != std::string::npos) {
            base_path = filename.substr(0, last_sep + 1);
        }
        
        DataImporter importer(base_path.empty() ? "./" : base_path);
        std::string file_only = (last_sep != std::string::npos) ? filename.substr(last_sep + 1) : filename;
        
        auto graph = importer.loadKAOGraph(file_only);
        if (!graph) {
            std::cerr << "Error: Could not load graph: " << filename << "\n";
            return 1;
        }

        // Auto-detect source and target from Targets array if -1
        if (s_vertex == -1 || t_vertex == -1) {
            int detected_s = -1, detected_t = -1;
            for (size_t i = 0; i < graph->target_vertices_.size(); ++i) {
                if (graph->target_vertices_[i] == 1) {
                    if (detected_s == -1) detected_s = static_cast<int>(i);
                    else if (detected_t == -1) detected_t = static_cast<int>(i);
                }
            }
            if (s_vertex == -1) s_vertex = detected_s;
            if (t_vertex == -1) t_vertex = detected_t;
            
            if (s_vertex == -1 || t_vertex == -1) {
                std::cerr << "Error: Could not auto-detect source/target vertices\n";
                return 1;
            }
            std::cout << "Auto-detected: s=" << s_vertex << ", t=" << t_vertex << "\n";
        }

        std::vector<std::string> method_names = {
            "Standard Factoring (Level 0)",
            "Recursive Decomposition (Level 1)", 
            "Simple Factoring (Level 2)",
            "M-Decomposition (Level 3)"
        };

        std::cout << "Running test on: " << filename << "\n";
        std::cout << "  Source: " << s_vertex << ", Target: " << t_vertex << "\n";
        std::cout << "  Diameter: " << diameter << ", Method: " << method_names[method_id] << "\n";
        std::cout << "  Repetitions: " << repetitions << "\n\n";

        double total_time = 0.0;
        double total_rel = 0.0;
        long long total_recs = 0;

        for (int rep = 0; rep < repetitions; ++rep) {
            // Reload graph for each repetition to ensure clean state
            auto g = importer.loadKAOGraph(file_only);
            
            ReliabilityResult result;
            switch (method_id) {
                case 0:  // Level 0: Global Standard Factoring (Baseline - SLOWEST)
                    result = g->calculateReliabilityBetweenVertices(s_vertex, t_vertex, diameter);
                    break;
                case 1:  // Level 1: Recursive Decomposition (Nested Recursion - INEFFICIENT)
                    result = g->calculateReliabilityWithRecursiveDecomposition(s_vertex, t_vertex, diameter);
                    break;
                case 2:  // Level 2: Iterative Decomp + Simple Facto (Convolution - FAST)
                    result = g->calculateReliabilityWithDecomposition(s_vertex, t_vertex, diameter);
                    break;
                case 3:  // Level 3: Iterative Decomp + Modified Facto (FASTEST)
                    result = g->calculateReliabilityWithMDecomposition(s_vertex, t_vertex, diameter);
                    break;
            }
            
            total_time += result.execution_time_sec;
            total_rel += result.reliability;
            total_recs += result.recursions;
            
            if (repetitions > 1) {
                std::cout << "  Rep " << (rep + 1) << ": Reliability=" << result.reliability 
                          << ", Time=" << result.execution_time_sec << "s\n";
            }
        }

        double avg_time = total_time / repetitions;
        double avg_rel = total_rel / repetitions;
        long long avg_recs = total_recs / repetitions;

        std::cout << "\nResults:\n";
        std::cout << "  Reliability: " << avg_rel << "\n";
        std::cout << "  Avg Time:    " << avg_time << " seconds\n";
        std::cout << "  Avg Recs:    " << avg_recs << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

/**
 * @brief Handle format conversion commands
 * @param args Command line arguments
 * @return Exit code
 */
int handleConversion(const std::vector<std::string>& args) {
    if (args.size() < 4) {
        std::cerr << "Error: Insufficient arguments for conversion\n";
        return 1;
    }

    const std::string& conversion_type = args[1];
    const std::string& input_file = args[2];
    const std::string& output_file = args[3];

    try {
        if (conversion_type == "edge2kao") {
            if (args.size() < 5) {
                std::cerr << "Error: Missing reliability value for edge2kao conversion\n";
                return 1;
            }
            double reliability = std::stod(args[4]);
            GraphOperations::convertEdgeListToKAOFO(input_file, output_file, reliability);
            std::cout << "Successfully converted Edge List to KAO format\n";
        } else if (conversion_type == "kao2edge") {
            GraphOperations::convertKAOFOToEdgeList(input_file, output_file);
            std::cout << "Successfully converted KAO format to Edge List\n";
        } else {
            std::cerr << "Error: Unknown conversion type: " << conversion_type << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during conversion: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

/**
 * @brief Handle test execution
 * @param args Command line arguments
 * @return Exit code
 */
int handleTesting(const std::vector<std::string>& args) {
    if (args.size() < 2) {
        std::cerr << "Error: Missing method ID for testing\n";
        return 1;
    }

    int method_id;
    try {
        method_id = std::stoi(args[1]);
        if (method_id < 0 || method_id > 3) {
            std::cerr << "Error: Method ID must be between 0 and 3\n";
            return 1;
        }
    } catch (const std::exception&) {
        std::cerr << "Error: Invalid method ID: " << args[1] << "\n";
        return 1;
    }

    std::string output_file = "test_results.csv";
    if (args.size() > 2) {
        output_file = args[2];
    }

    try {
        DataImporter importer("graphs_data/");
        TestSuite test_suite(importer);

        // Set up test configurations
        std::map<std::string, std::vector<int>> test_configs = {
            {"3_blocks_sausage_3x3_kao.txt", {9, 10, 11, 12}},
            {"4_blocks_sausage_3x3_kao.txt", {13, 14, 15, 16}},
            {"5_blocks_sausage_3x3_kao.txt", {17, 18, 19, 20}},
            {"6_blocks_sausage_3x3_kao.txt", {21, 22, 23, 24}}
        };
        test_suite.setTestConfigurations(test_configs);

        std::cout << "Running comprehensive tests with method " << method_id << "...\n";
        bool success = test_suite.runComprehensiveTests(method_id, output_file);
        
        if (success) {
            std::cout << "All tests completed successfully. Results saved to " << output_file << "\n";
        } else {
            std::cerr << "Some tests failed. Check the output for details.\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during testing: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

/**
 * @brief Main application entry point
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code
 */
int main(int argc, char* argv[]) {
    std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty() || args[0] == "--help") {
        printUsage();
        return 0;
    }

    try {
        if (args[0] == "--convert") {
            return handleConversion(args);
        } else if (args[0] == "--test") {
            return handleTesting(args);
        } else if (args[0] == "--run") {
            return handleRun(args);
        } else {
            std::cerr << "Error: Unknown command: " << args[0] << "\n";
            std::cerr << "Use --help for usage information.\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
