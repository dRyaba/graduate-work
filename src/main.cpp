/**
 * @file main.cpp
 * @brief Main application for Graph Reliability Analysis
 * @author Graduate Work Project
 * @date 2026
 * 
 * This is the main entry point for the graph reliability analysis application.
 * It provides command-line interface for running reliability tests and
 * graph format conversions.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include "graph_reliability.h"
#include "graph_reliability/Logger.h"
#include "graph_reliability/GraphVisualizer.h"

using namespace graph_reliability;

/**
 * @brief Print application usage information
 */
void printUsage() {
    std::cout << "Graph Reliability Analysis Tool v" << getVersion() << "\n\n";
    std::cout << "Usage:\n";
#ifdef _WIN32
    std::cout << "  .\\build\\graph_reliability.exe [options]\n\n";
#else
    std::cout << "  ./build/graph_reliability [options]\n\n";
#endif
    std::cout << "Options:\n";
    std::cout << "  --run <file> <s> <t> <d> <method> [reps]           Run test on specific graph file\n";
    std::cout << "  --test <method_id> [output_file]                   Run comprehensive tests on predefined graphs\n";
    std::cout << "  --convert edge2kao <input> <output> <reliability>  Convert Edge List to KAO format\n";
    std::cout << "  --convert kao2edge <input> <output>                Convert KAO format to Edge List\n";
    std::cout << "  --visualize <graph_file> <output_file> [opts]      Visualize graph (SVG or DOT)\n";
    std::cout << "  --help                                             Show this help message\n\n";
    std::cout << "Test Methods (aligned with Optimization Levels):\n";
    std::cout << "  0 - Standard Factoring      (Level 0: Global, No Decomposition - SLOWEST)\n";
    std::cout << "  1 - Recursive Decomposition (Level 1: Nested Recursion - INEFFICIENT)\n";
    std::cout << "  2 - Simple Factoring        (Level 2: Convolution + Simple Facto - FAST)\n";
    std::cout << "  3 - M-Decomposition         (Level 3: Convolution + Modified Facto - FASTEST)\n";
    std::cout << "  4 - Cancela-Petingi         (Level 4: Path-based factoring with SPT)\n";
    std::cout << "  5 - M-Decomp + CPFM         (Level 5: Decomposition + Path-based factoring)\n\n";
    std::cout << "Parameters for --run:\n";
    std::cout << "  <file>    - Graph file path (KAO format)\n";
    std::cout << "  <s>       - Source vertex (0-based index, or -1 to auto-detect from targets)\n";
    std::cout << "  <t>       - Target vertex (0-based index, or -1 to auto-detect from targets)\n";
    std::cout << "  --1-based - If present before <file>, s and t are interpreted as 1-based (as in VKR/course work)\n";
    std::cout << "  <d>       - Diameter upper bound\n";
    std::cout << "  <method>  - Method ID (0-5)\n";
    std::cout << "  [reps]    - Number of repetitions (default: 1)\n\n";
    std::cout << "Examples:\n";
#ifdef _WIN32
    std::cout << "  .\\build\\graph_reliability.exe --run graphs_data\\K4_kao.txt 0 3 3 4 1\n";
    std::cout << "  .\\build\\graph_reliability.exe --run graphs_data\\3_blocks_sausage_3x3_kao.txt -1 -1 12 4 1\n";
    std::cout << "  .\\build\\graph_reliability.exe --test 3 results.csv\n";
    std::cout << "  .\\build\\graph_reliability.exe --convert edge2kao input.edgelist output.kao 0.9\n";
    std::cout << "  .\\build\\graph_reliability.exe --visualize graphs_data\\K4_kao.txt k4.svg --probs\n";
    std::cout << "  .\\build\\graph_reliability.exe --visualize graphs_data\\Geant2004_kao.txt geant.svg --s 11 --t 95\n";
    std::cout << "  .\\build\\graph_reliability.exe --visualize graphs_data\\Geant2004_kao.txt geant.dot --dot\n";
#else
    std::cout << "  ./build/graph_reliability --run graphs_data/K4_kao.txt 0 3 3 4 1\n";
    std::cout << "  ./build/graph_reliability --run graphs_data/3_blocks_sausage_3x3_kao.txt -1 -1 12 4 1\n";
    std::cout << "  ./build/graph_reliability --test 3 results.csv\n";
    std::cout << "  ./build/graph_reliability --convert edge2kao input.edgelist output.kao 0.9\n";
    std::cout << "  ./build/graph_reliability --visualize graphs_data/K4_kao.txt k4.svg --probs\n";
    std::cout << "  ./build/graph_reliability --visualize graphs_data/Geant2004_kao.txt geant.svg --s 11 --t 95\n";
    std::cout << "  ./build/graph_reliability --visualize graphs_data/Geant2004_kao.txt geant.dot --dot\n";
#endif
    std::cout << "\nParameters for --visualize:\n";
    std::cout << "  <graph_file>  - Input graph (KAO format)\n";
    std::cout << "  <output_file> - Output path (.svg or .dot)\n";
    std::cout << "  --s <n>       - Source vertex index (highlighted green)\n";
    std::cout << "  --t <n>       - Target vertex index (highlighted red)\n";
    std::cout << "  --dot         - Write Graphviz DOT instead of SVG\n";
    std::cout << "  --probs       - Label edges with probabilities\n";
    std::cout << "  --iter <n>    - FR layout iterations (default 500)\n";
}

/**
 * @brief Handle single graph test run
 * @param args Command line arguments
 * @return Exit code
 */
int handleRun(const std::vector<std::string>& args) {
    // --run [--1-based] <file> <s> <t> <d> <method> [reps]
    bool one_based = false;
    size_t arg_offset = 1;
    if (args.size() >= 7 && args[1] == "--1-based") {
        one_based = true;
        arg_offset = 2;
    }
    if (args.size() < 5 + arg_offset) {
        LOG_ERROR("Insufficient arguments for --run");
        std::cerr << "Error: Insufficient arguments for --run\n";
        std::cerr << "Usage: --run [--1-based] <file> <s> <t> <d> <method> [reps]\n";
        return 1;
    }

    const std::string& filename = args[arg_offset];
    int s_vertex, t_vertex, diameter, method_id;
    int repetitions = 1;

    try {
        s_vertex = std::stoi(args[arg_offset + 1]);
        t_vertex = std::stoi(args[arg_offset + 2]);
        diameter = std::stoi(args[arg_offset + 3]);
        method_id = std::stoi(args[arg_offset + 4]);
        
        if (args.size() > arg_offset + 5) {
            repetitions = std::stoi(args[arg_offset + 5]);
        }
        
        if (one_based && s_vertex != -1 && t_vertex != -1) {
            s_vertex--;
            t_vertex--;
        }

        if (method_id < 0 || method_id > 5) {
            LOG_ERROR("Method ID must be between 0 and 5, got: {}", method_id);
            std::cerr << "Error: Method ID must be between 0 and 5\n";
            return 1;
        }
        if (diameter < 1) {
            LOG_ERROR("Diameter must be positive, got: {}", diameter);
            std::cerr << "Error: Diameter must be positive\n";
            return 1;
        }
        if (repetitions < 1) {
            LOG_ERROR("Repetitions must be positive, got: {}", repetitions);
            std::cerr << "Error: Repetitions must be positive\n";
            return 1;
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Invalid numeric argument: {}", e.what());
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
        
        LOG_INFO("Loading graph from: {}", filename);
        auto graph = importer.loadKAOGraph(file_only);
        if (!graph) {
            LOG_ERROR("Could not load graph: {}", filename);
            std::cerr << "Error: Could not load graph: " << filename << "\n";
            return 1;
        }
        LOG_DEBUG("Graph loaded successfully: {} vertices, {} edges", 
                  graph->numVertices(), graph->numEdges());

        int num_vertices = static_cast<int>(graph->numVertices());
        
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
                LOG_ERROR("Could not auto-detect source/target vertices");
                std::cerr << "Error: Could not auto-detect source/target vertices\n";
                return 1;
            }
            LOG_INFO("Auto-detected source={}, target={}", s_vertex, t_vertex);
            std::cout << "Auto-detected: s=" << s_vertex << ", t=" << t_vertex << "\n";
        }
        
        // Validate vertex indices
        if (s_vertex < 0 || s_vertex >= num_vertices) {
            std::cerr << "Error: Source vertex " << s_vertex << " is out of range [0, " 
                      << (num_vertices - 1) << "]\n";
            return 1;
        }
        if (t_vertex < 0 || t_vertex >= num_vertices) {
            std::cerr << "Error: Target vertex " << t_vertex << " is out of range [0, " 
                      << (num_vertices - 1) << "]\n";
            return 1;
        }
        if (s_vertex == t_vertex) {
            std::cerr << "Error: Source and target vertices must be different\n";
            return 1;
        }

        std::vector<std::string> method_names = {
            "Standard Factoring (Level 0)",
            "Recursive Decomposition (Level 1)", 
            "Simple Factoring (Level 2)",
            "M-Decomposition (Level 3)",
            "Cancela-Petingi (Level 4)",
            "M-Decomp + CPFM (Level 5)"
        };

        LOG_INFO("Starting reliability calculation: file={}, s={}, t={}, diameter={}, method={}, reps={}",
                 filename, s_vertex, t_vertex, diameter, method_names[method_id], repetitions);
        std::cout << "Running test on: " << filename << "\n";
        std::cout << "  Source: " << s_vertex << ", Target: " << t_vertex << "\n";
        std::cout << "  Diameter: " << diameter << ", Method: " << method_names[method_id] << "\n";
        std::cout << "  Repetitions: " << repetitions << "\n\n";

        double total_time = 0.0;
        double total_rel = 0.0;
        long long total_recs = 0;

        for (int rep = 0; rep < repetitions; ++rep) {
            LOG_DEBUG("Repetition {}/{}", rep + 1, repetitions);
            // Reload graph for each repetition to ensure clean state
            auto g = importer.loadKAOGraph(file_only);
            
            ReliabilityResult result;
            switch (method_id) {
                case 0:  // Level 0: Global Standard Factoring (Baseline - SLOWEST)
                    LOG_DEBUG("Using Standard Factoring method");
                    result = g->calculateReliabilityBetweenVertices(s_vertex, t_vertex, diameter);
                    break;
                case 1:  // Level 1: Recursive Decomposition (Nested Recursion - INEFFICIENT)
                    LOG_DEBUG("Using Recursive Decomposition method");
                    result = g->calculateReliabilityWithRecursiveDecomposition(s_vertex, t_vertex, diameter);
                    break;
                case 2:  // Level 2: Iterative Decomp + Simple Facto (Convolution - FAST)
                    LOG_DEBUG("Using Simple Factoring method");
                    result = g->calculateReliabilityWithDecomposition(s_vertex, t_vertex, diameter);
                    break;
                case 3:  // Level 3: Iterative Decomp + Modified Facto (FASTEST)
                    LOG_DEBUG("Using M-Decomposition method");
                    result = g->calculateReliabilityWithMDecomposition(s_vertex, t_vertex, diameter);
                    break;
                case 4:  // Level 4: Cancela-Petingi path-based factoring
                    LOG_DEBUG("Using Cancela-Petingi method");
                    result = g->calculateReliabilityCancelaPetingi(s_vertex, t_vertex, diameter);
                    break;
                case 5:  // Level 5: M-Decomposition + Cancela-Petingi hybrid
                    LOG_DEBUG("Using M-Decomp + CPFM method");
                    result = g->calculateReliabilityWithMDecompositionCPFM(s_vertex, t_vertex, diameter);
                    break;
            }
            
            LOG_INFO("Repetition {} completed: reliability={}, time={}s, recursions={}",
                     rep + 1, result.reliability, result.execution_time_sec, result.recursions);
            
            total_time += result.execution_time_sec;
            total_rel += result.reliability;
            total_recs += result.recursions;
            
            if (repetitions > 1) {
                std::cout << "  Rep " << (rep + 1) << ": Reliability=" << std::fixed << std::setprecision(15) << result.reliability 
                          << ", Time=" << result.execution_time_sec << "s\n";
            }
        }

        double avg_time = total_time / repetitions;
        double avg_rel = total_rel / repetitions;
        long long avg_recs = total_recs / repetitions;

        LOG_INFO("All repetitions completed: avg_reliability={}, avg_time={}s, avg_recursions={}",
                  avg_rel, avg_time, avg_recs);
        std::cout << "\nResults:\n";
        std::cout << "  Reliability: " << std::fixed << std::setprecision(15) << avg_rel << "\n";
        std::cout << "  Avg Time:    " << avg_time << " seconds\n";
        std::cout << "  Avg Recs:    " << avg_recs << "\n";

    } catch (const std::exception& e) {
        LOG_ERROR("Exception in handleRun: {}", e.what());
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
        LOG_ERROR("Insufficient arguments for conversion");
        std::cerr << "Error: Insufficient arguments for conversion\n";
        return 1;
    }

    const std::string& conversion_type = args[1];
    const std::string& input_file = args[2];
    const std::string& output_file = args[3];

    LOG_INFO("Starting conversion: type={}, input={}, output={}", conversion_type, input_file, output_file);

    try {
        if (conversion_type == "edge2kao") {
            if (args.size() < 5) {
                LOG_ERROR("Missing reliability value for edge2kao conversion");
                std::cerr << "Error: Missing reliability value for edge2kao conversion\n";
                return 1;
            }
            double reliability = std::stod(args[4]);
            LOG_DEBUG("Converting Edge List to KAO with reliability={}", reliability);
            GraphOperations::convertEdgeListToKAOFO(input_file, output_file, reliability);
            LOG_INFO("Successfully converted Edge List to KAO format");
            std::cout << "Successfully converted Edge List to KAO format\n";
        } else if (conversion_type == "kao2edge") {
            LOG_DEBUG("Converting KAO to Edge List");
            GraphOperations::convertKAOFOToEdgeList(input_file, output_file);
            LOG_INFO("Successfully converted KAO format to Edge List");
            std::cout << "Successfully converted KAO format to Edge List\n";
        } else {
            LOG_ERROR("Unknown conversion type: {}", conversion_type);
            std::cerr << "Error: Unknown conversion type: " << conversion_type << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Exception during conversion: {}", e.what());
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
        if (method_id < 0 || method_id > 5) {
            std::cerr << "Error: Method ID must be between 0 and 5\n";
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

        LOG_INFO("Running comprehensive tests with method {}", method_id);
        std::cout << "Running comprehensive tests with method " << method_id << "...\n";
        bool success = test_suite.runComprehensiveTests(method_id, output_file);
        
        if (success) {
            LOG_INFO("All tests completed successfully. Results saved to {}", output_file);
            std::cout << "All tests completed successfully. Results saved to " << output_file << "\n";
        } else {
            LOG_ERROR("Some tests failed");
            std::cerr << "Some tests failed. Check the output for details.\n";
            return 1;
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Exception during testing: {}", e.what());
        std::cerr << "Error during testing: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

/**
 * @brief Handle --visualize command: layout + SVG/DOT export
 * @param args Command line arguments (args[0] == "--visualize")
 * @return Exit code
 */
int handleVisualize(const std::vector<std::string>& args) {
    // --visualize <graph_file> <output_file> [--s <n>] [--t <n>] [--dot] [--probs] [--iter <n>]
    if (args.size() < 3) {
        std::cerr << "Error: Insufficient arguments for --visualize\n";
        std::cerr << "Usage: --visualize <graph_file> <output_file> [--s <n>] [--t <n>] [--dot] [--probs] [--iter <n>]\n";
        return 1;
    }

    const std::string& filename    = args[1];
    const std::string& output_path = args[2];

    GraphVisualizer::Options opts;
    bool output_dot = false;

    for (size_t i = 3; i < args.size(); ++i) {
        if (args[i] == "--dot") {
            output_dot = true;
        } else if (args[i] == "--probs") {
            opts.show_probs = true;
        } else if (args[i] == "--s" && i + 1 < args.size()) {
            opts.source = std::stoi(args[++i]);
        } else if (args[i] == "--t" && i + 1 < args.size()) {
            opts.target = std::stoi(args[++i]);
        } else if (args[i] == "--iter" && i + 1 < args.size()) {
            opts.iterations = std::stoi(args[++i]);
        }
    }

    try {
        // Load graph
        std::string base_path;
        size_t last_sep = filename.find_last_of("/\\");
        if (last_sep != std::string::npos)
            base_path = filename.substr(0, last_sep + 1);

        DataImporter importer(base_path.empty() ? "./" : base_path);
        std::string file_only = (last_sep != std::string::npos)
                                ? filename.substr(last_sep + 1) : filename;

        auto graph = importer.loadKAOGraph(file_only);
        if (!graph) {
            std::cerr << "Error: Could not load graph: " << filename << "\n";
            return 1;
        }

        // Auto-detect s/t from target_vertices_ if not provided
        if (opts.source < 0 || opts.target < 0) {
            int detected_s = -1, detected_t = -1;
            for (size_t i = 0; i < graph->target_vertices_.size(); ++i) {
                if (graph->target_vertices_[i] == 1) {
                    if (detected_s < 0) detected_s = static_cast<int>(i);
                    else if (detected_t < 0) detected_t = static_cast<int>(i);
                }
            }
            if (opts.source < 0) opts.source = detected_s;
            if (opts.target < 0) opts.target = detected_t;
        }

        if (output_dot)
            GraphVisualizer::exportDot(*graph, output_path, opts);
        else
            GraphVisualizer::exportSVG(*graph, output_path, opts);

        std::cout << "Exported " << graph->numVertices() << " vertices, "
                  << (graph->fo_.size() / 2) << " edges"
                  << " → " << output_path << "\n";
        if (opts.source >= 0)
            std::cout << "  s=" << opts.source << " (green)  t=" << opts.target << " (red)\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
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
    // Initialize logger if logging is enabled
    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    // Check for verbose flag
    bool verbose = false;
    std::vector<std::string> args(argv + 1, argv + argc);
    for (const auto& arg : args) {
        if (arg == "--verbose" || arg == "-v") {
            verbose = true;
            break;
        }
    }
    
    graph_reliability::Logger::initialize();
    if (verbose) {
        graph_reliability::Logger::setLevel("debug");
    } else {
        graph_reliability::Logger::setLevel("info");
    }
    LOG_INFO("Graph Reliability Analysis Tool v{} started", graph_reliability::getVersion());
    #else
    std::vector<std::string> args(argv + 1, argv + argc);
    #endif

    if (args.empty() || args[0] == "--help") {
        printUsage();
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        graph_reliability::Logger::shutdown();
        #endif
        return 0;
    }

    try {
        int result = 0;
        if (args[0] == "--convert") {
            result = handleConversion(args);
        } else if (args[0] == "--test") {
            result = handleTesting(args);
        } else if (args[0] == "--run") {
            result = handleRun(args);
        } else if (args[0] == "--visualize") {
            result = handleVisualize(args);
        } else {
            LOG_ERROR("Unknown command: {}", args[0]);
            std::cerr << "Error: Unknown command: " << args[0] << "\n";
            std::cerr << "Use --help for usage information.\n";
            result = 1;
        }
        
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        LOG_INFO("Application exiting with code {}", result);
        graph_reliability::Logger::shutdown();
        #endif
        
        return result;
    } catch (const std::exception& e) {
        LOG_CRITICAL("Fatal error: {}", e.what());
        std::cerr << "Fatal error: " << e.what() << "\n";
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        graph_reliability::Logger::shutdown();
        #endif
        return 2;
    }
}
