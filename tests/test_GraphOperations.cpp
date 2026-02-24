/**
 * @file test_GraphOperations.cpp
 * @brief Unit tests for GraphOperations class
 * @author Graduate Work Project
 * @date 2026
 */

#include <gtest/gtest.h>
#include "graph_reliability/GraphOperations.h"
#include "graph_reliability/DataImporter.h"
#include "graph_reliability/Exceptions.h"
#include <filesystem>

using namespace graph_reliability;

class GraphOperationsTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_data_dir_ = "test_data_temp";
        std::filesystem::create_directories(test_data_dir_);
    }
    
    void TearDown() override {
        if (std::filesystem::exists(test_data_dir_)) {
            std::filesystem::remove_all(test_data_dir_);
        }
    }
    
    std::string test_data_dir_;
};

TEST_F(GraphOperationsTest, ConvertEdgeListToKAO) {
    // Create a simple edge list file
    std::string input_file = test_data_dir_ + "/input.edgelist";
    std::ofstream file(input_file);
    file << "1 -- 2\n";
    file << "2 -- 3\n";
    file << "3 -- 1\n";
    file.close();
    
    std::string output_file = test_data_dir_ + "/output.kao";
    
    GraphOperations::convertEdgeListToKAOFO(input_file, output_file, 0.9);
    
    // Check that output file was created
    EXPECT_TRUE(std::filesystem::exists(output_file));
    
    // Verify file is not empty
    std::ifstream check(output_file);
    std::string line;
    std::getline(check, line);
    EXPECT_FALSE(line.empty());
}

TEST_F(GraphOperationsTest, ConvertKAOToEdgeList) {
    // Create a simple KAO file
    std::string input_file = test_data_dir_ + "/input.kao";
    std::ofstream file(input_file);
    file << "0,2,4,6\n"; // KAO
    file << "1,2,0,2,0,1\n"; // FO (1-based)
    file << "0,1,0,1\n"; // Targets
    file << "0.9\n"; // Probability
    file.close();
    
    std::string output_file = test_data_dir_ + "/output.edgelist";
    
    GraphOperations::convertKAOFOToEdgeList(input_file, output_file);
    
    // Check that output file was created
    EXPECT_TRUE(std::filesystem::exists(output_file));
    
    // Verify file contains edges
    std::ifstream check(output_file);
    std::string line;
    bool has_content = false;
    while (std::getline(check, line)) {
        if (!line.empty() && line.find("--") != std::string::npos) {
            has_content = true;
            break;
        }
    }
    EXPECT_TRUE(has_content);
}

TEST_F(GraphOperationsTest, ConvertThrowsOnMissingFile) {
    EXPECT_THROW(
        GraphOperations::convertEdgeListToKAOFO("nonexistent.edgelist", "output.kao", 0.9),
        FileNotFoundException
    );
    
    EXPECT_THROW(
        GraphOperations::convertKAOFOToEdgeList("nonexistent.kao", "output.edgelist"),
        FileNotFoundException
    );
}
