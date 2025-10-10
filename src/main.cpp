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
    std::cout << "  --convert edge2kao <input> <output> <reliability>  Convert Edge List to KAO format\n";
    std::cout << "  --convert kao2edge <input> <output>                Convert KAO format to Edge List\n";
    std::cout << "  --test <method_id> [output_file]                   Run comprehensive tests\n";
    std::cout << "  --help                                            Show this help message\n\n";
    std::cout << "Test Methods:\n";
    std::cout << "  0 - M-Decomposition\n";
    std::cout << "  1 - Recursive Decomposition\n";
    std::cout << "  2 - Simple Factoring\n";
    std::cout << "  3 - Standard Factoring\n\n";
    std::cout << "Examples:\n";
    std::cout << "  ./graph_reliability --convert edge2kao input.edgelist output.kao 0.9\n";
    std::cout << "  ./graph_reliability --convert kao2edge input.kao output.edgelist\n";
    std::cout << "  ./graph_reliability --test 0 results.csv\n";
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
