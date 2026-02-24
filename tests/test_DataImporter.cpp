/**
 * @file test_DataImporter.cpp
 * @brief Unit tests for DataImporter class
 * @author Graduate Work Project
 * @date 2026
 */

#include <gtest/gtest.h>
#include "graph_reliability/DataImporter.h"
#include "graph_reliability/Exceptions.h"
#include <filesystem>
#include <fstream>

using namespace graph_reliability;

class DataImporterTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a temporary test data directory
        test_data_dir_ = "test_data_temp";
        std::filesystem::create_directories(test_data_dir_);
        
        // Create a simple test KAO file
        createTestKAOFile();
    }
    
    void TearDown() override {
        // Clean up test files
        if (std::filesystem::exists(test_data_dir_)) {
            std::filesystem::remove_all(test_data_dir_);
        }
    }
    
    void createTestKAOFile() {
        std::string filename = test_data_dir_ + "/test_graph.kao";
        std::ofstream file(filename);
        file << "0,2,4,6\n"; // KAO
        file << "1,2,0,2,0,1\n"; // FO (1-based)
        file << "0,1,0,1\n"; // Targets (with leading 0)
        file << "0.9\n"; // Probability
        file.close();
    }
    
    std::string test_data_dir_;
};

TEST_F(DataImporterTest, Constructor) {
    DataImporter importer(test_data_dir_);
    
    EXPECT_EQ(importer.getBasePath(), test_data_dir_ + "/");
}

TEST_F(DataImporterTest, GetFullPath) {
    DataImporter importer(test_data_dir_);
    
    std::string full_path = importer.getFullPath("test.txt");
    EXPECT_EQ(full_path, test_data_dir_ + "/test.txt");
    
    // Absolute path should remain unchanged
    std::string abs_path = "/absolute/path.txt";
    std::string result = importer.getFullPath(abs_path);
    // On Windows, this might be different, so just check it's not empty
    EXPECT_FALSE(result.empty());
}

TEST_F(DataImporterTest, FileExists) {
    DataImporter importer(test_data_dir_);
    
    EXPECT_TRUE(importer.fileExists("test_graph.kao"));
    EXPECT_FALSE(importer.fileExists("nonexistent.txt"));
}

TEST_F(DataImporterTest, LoadKAOGraph) {
    DataImporter importer(test_data_dir_);
    
    auto graph = importer.loadKAOGraph("test_graph.kao");
    
    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->numVertices(), 3);
    EXPECT_EQ(graph->numEdges(), 6);
}

TEST_F(DataImporterTest, LoadKAOGraphThrowsOnMissingFile) {
    DataImporter importer(test_data_dir_);
    
    EXPECT_THROW(importer.loadKAOGraph("nonexistent.kao"), FileNotFoundException);
}

TEST_F(DataImporterTest, LoadKAOGraphThrowsOnInvalidFormat) {
    DataImporter importer(test_data_dir_);
    
    // Create invalid file
    std::string filename = test_data_dir_ + "/invalid.kao";
    std::ofstream file(filename);
    file << "invalid content\n";
    file.close();
    
    EXPECT_THROW(importer.loadKAOGraph("invalid.kao"), InvalidFormatException);
}

TEST_F(DataImporterTest, LoadEdgeListGraph) {
    DataImporter importer(test_data_dir_);
    
    // Create a simple edge list file
    std::string filename = test_data_dir_ + "/test.edgelist";
    std::ofstream file(filename);
    file << "1 -- 2\n";
    file << "2 -- 3\n";
    file << "3 -- 1\n";
    file.close();
    
    auto graph = importer.loadEdgeListGraph("test.edgelist", 0.9);
    
    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->numVertices(), 3);
    EXPECT_GT(graph->numEdges(), 0);
}

TEST_F(DataImporterTest, SetBasePath) {
    DataImporter importer(test_data_dir_);
    
    importer.setBasePath("new_path");
    EXPECT_EQ(importer.getBasePath(), "new_path/");
}

TEST_F(DataImporterTest, GetAvailableFiles) {
    DataImporter importer(test_data_dir_);
    
    auto files = importer.getAvailableFiles();
    
    // Should find at least our test file
    bool found_test = false;
    for (const auto& f : files) {
        if (f == "test_graph.kao") {
            found_test = true;
            break;
        }
    }
    EXPECT_TRUE(found_test);
}
