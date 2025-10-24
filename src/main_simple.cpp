#include <iostream>
#include <string>
#include <vector>
#include "DataImporter.h"
#include "kGraphOperations.h"

// Forward declarations for conversion functions
void convertEdgeListToKAOFO(const std::string& inPath, const std::string& outPath, double reliability);
void convertKAOFOToEdgeList(const std::string& inPath, const std::string& outPath);

int main(int argc, char* argv[]) {
    std::cout << "Graph Reliability Analysis Tool v1.0.0\n\n";
    
    if (argc < 2) {
        std::cout << "Usage:\n";
        std::cout << "  ./graph_reliability --help                    Show this help message\n";
        std::cout << "  ./graph_reliability --test <method_id>        Run tests\n";
        std::cout << "  ./graph_reliability --convert <type> <args>   Convert formats\n\n";
        std::cout << "Test Methods:\n";
        std::cout << "  0 - M-Decomposition\n";
        std::cout << "  1 - Recursive Decomposition\n";
        std::cout << "  2 - Simple Factoring\n";
        std::cout << "  3 - Standard Factoring\n\n";
        std::cout << "Conversion Examples:\n";
        std::cout << "  ./graph_reliability --convert edge2kao input.edgelist output.kao 0.9\n";
        std::cout << "  ./graph_reliability --convert kao2edge input.kao output.edgelist\n";
        return 0;
    }
    
    std::string command = argv[1];
    
    if (command == "--help") {
        std::cout << "Graph Reliability Analysis Tool v1.0.0\n\n";
        std::cout << "Usage:\n";
        std::cout << "  ./graph_reliability --help                    Show this help message\n";
        std::cout << "  ./graph_reliability --test <method_id>        Run tests\n";
        std::cout << "  ./graph_reliability --convert <type> <args>   Convert formats\n\n";
        std::cout << "Test Methods:\n";
        std::cout << "  0 - M-Decomposition\n";
        std::cout << "  1 - Recursive Decomposition\n";
        std::cout << "  2 - Simple Factoring\n";
        std::cout << "  3 - Standard Factoring\n\n";
        std::cout << "Conversion Examples:\n";
        std::cout << "  ./graph_reliability --convert edge2kao input.edgelist output.kao 0.9\n";
        std::cout << "  ./graph_reliability --convert kao2edge input.kao output.edgelist\n";
        return 0;
    }
    
    if (command == "--test") {
        if (argc < 3) {
            std::cerr << "Error: Missing method ID for testing\n";
            return 1;
        }
        
        int method_id = std::stoi(argv[2]);
        std::cout << "Running tests with method " << method_id << "...\n";
        
        try {
            DataImporter importer("graphs_data/");
            
            // Test with a simple graph
            auto graph = importer.loadKAOGraph("K4_kao.txt");
            if (graph) {
                std::cout << "Successfully loaded K4 graph\n";
                std::cout << "Graph has " << graph->KAO.size() << " vertices\n";
            } else {
                std::cerr << "Failed to load test graph\n";
                return 1;
            }
            
            std::cout << "Test completed successfully!\n";
        } catch (const std::exception& e) {
            std::cerr << "Error during testing: " << e.what() << "\n";
            return 1;
        }
        
        return 0;
    }
    
    if (command == "--convert") {
        if (argc < 4) {
            std::cerr << "Error: Insufficient arguments for conversion\n";
            return 1;
        }
        
        std::string conversion_type = argv[2];
        std::string input_file = argv[3];
        std::string output_file = argv[4];
        
        try {
            if (conversion_type == "edge2kao") {
                if (argc < 6) {
                    std::cerr << "Error: Missing reliability value for edge2kao conversion\n";
                    return 1;
                }
                double reliability = std::stod(argv[5]);
                // Use the existing conversion function
                convertEdgeListToKAOFO(input_file, output_file, reliability);
                std::cout << "Successfully converted Edge List to KAO format\n";
            } else if (conversion_type == "kao2edge") {
                // Use the existing conversion function
                convertKAOFOToEdgeList(input_file, output_file);
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
    
    std::cerr << "Error: Unknown command: " << command << "\n";
    std::cerr << "Use --help for usage information.\n";
    return 1;
}
