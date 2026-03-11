/**
 * @file test_ReliabilityGraph.cpp
 * @brief Unit tests for ReliabilityGraph class
 * @author Graduate Work Project
 * @date 2026
 */

#include <gtest/gtest.h>
#include "graph_reliability/ReliabilityGraph.h"
#include "graph_reliability/DataImporter.h"
#include "graph_reliability/Exceptions.h"

using namespace graph_reliability;

class ReliabilityGraphTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple graph: triangle (3 vertices, 3 edges)
        kao_ = {0, 2, 4, 6};
        fo_ = {1, 2, 0, 2, 0, 1};
        p_array_ = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
        targets_ = {1, 0, 1}; // Vertices 0 and 2 are targets
    }

    std::vector<ReliabilityGraph::VertexId> kao_;
    std::vector<ReliabilityGraph::VertexId> fo_;
    std::vector<ReliabilityGraph::Probability> p_array_;
    std::vector<int> targets_;
};

TEST_F(ReliabilityGraphTest, Constructor) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    EXPECT_EQ(graph.numVertices(), 3);
    EXPECT_EQ(graph.numEdges(), 6);
    EXPECT_EQ(graph.target_vertices_.size(), 3);
    EXPECT_EQ(graph.target_vertices_[0], 1);
    EXPECT_EQ(graph.target_vertices_[2], 1);
}

TEST_F(ReliabilityGraphTest, IsKConnected) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    // Triangle with targets at vertices 0 and 2 should be connected
    EXPECT_TRUE(graph.isKConnected());
    
    // Single target
    std::vector<int> single_target = {1, 0, 0};
    ReliabilityGraph graph2(kao_, fo_, p_array_, single_target);
    EXPECT_TRUE(graph2.isKConnected());
}

TEST_F(ReliabilityGraphTest, FindLastUnreliableEdge) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    auto edge = graph.findLastUnreliableEdge();
    EXPECT_TRUE(edge.has_value());
    EXPECT_LT(*edge, graph.numEdges());
    
    // All edges reliable
    std::vector<ReliabilityGraph::Probability> all_reliable(6, 1.0);
    ReliabilityGraph reliable_graph(kao_, fo_, all_reliable, targets_);
    auto no_edge = reliable_graph.findLastUnreliableEdge();
    EXPECT_FALSE(no_edge.has_value());
}

TEST_F(ReliabilityGraphTest, RemoveEdge) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    ReliabilityGraph modified = graph.removeEdge(0, 1);
    
    EXPECT_FALSE(modified.hasEdge(0, 1));
    EXPECT_FALSE(modified.hasEdge(1, 0));
    EXPECT_EQ(modified.target_vertices_.size(), graph.target_vertices_.size());
}

TEST_F(ReliabilityGraphTest, CalculateReliabilityBetweenVertices) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    // Calculate reliability between vertices 0 and 2 with diameter 2
    auto result = graph.calculateReliabilityBetweenVertices(0, 2, 2);
    
    EXPECT_GE(result.reliability, 0.0);
    EXPECT_LE(result.reliability, 1.0);
    EXPECT_GT(result.recursions, 0);
    EXPECT_GT(result.execution_time_sec, 0.0);
}

TEST_F(ReliabilityGraphTest, CalculateReliabilityWithMDecomposition) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    // For a simple triangle, M-decomposition should work
    auto result = graph.calculateReliabilityWithMDecomposition(0, 2, 2);
    
    EXPECT_GE(result.reliability, 0.0);
    EXPECT_LE(result.reliability, 1.0);
    EXPECT_GE(result.recursions, 0);
    EXPECT_GT(result.execution_time_sec, 0.0);
}

TEST_F(ReliabilityGraphTest, CalculateReliabilityWithDecomposition) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    auto result = graph.calculateReliabilityWithDecomposition(0, 2, 2);
    
    EXPECT_GE(result.reliability, 0.0);
    EXPECT_LE(result.reliability, 1.0);
    EXPECT_GE(result.recursions, 0);
    EXPECT_GT(result.execution_time_sec, 0.0);
}

TEST_F(ReliabilityGraphTest, CalculateReliabilityCancelaPetingi) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    auto result = graph.calculateReliabilityCancelaPetingi(0, 2, 2);
    
    EXPECT_GE(result.reliability, 0.0);
    EXPECT_LE(result.reliability, 1.0);
    EXPECT_GE(result.recursions, 0);
    EXPECT_GT(result.execution_time_sec, 0.0);
}

TEST_F(ReliabilityGraphTest, CancelaPetingiMatchesStandardFactoring) {
    ReliabilityGraph graph(kao_, fo_, p_array_, targets_);
    
    auto cp_result = graph.calculateReliabilityCancelaPetingi(0, 2, 2);
    auto std_result = graph.calculateReliabilityBetweenVertices(0, 2, 2);
    
    EXPECT_NEAR(cp_result.reliability, std_result.reliability, 1e-10);
}

TEST_F(ReliabilityGraphTest, CancelaPetingiK4RecursionCount) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("K4_kao.txt")) {
        GTEST_SKIP() << "K4_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("K4_kao.txt");
    ASSERT_TRUE(graph);

    auto result = graph->calculateReliabilityCancelaPetingi(0, 3, 3);

    EXPECT_NEAR(result.reliability, 0.997848, 1e-10);
    EXPECT_EQ(result.recursions, 17);  // Optimized with corrected ISPT
}

TEST_F(ReliabilityGraphTest, CancelaPetingiK5RecursionCount) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("K5_kao.txt")) {
        GTEST_SKIP() << "K5_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("K5_kao.txt");
    ASSERT_TRUE(graph);

    auto result1 = graph->calculateReliabilityCancelaPetingi(0, 4, 4);
    auto result2 = graph->calculateReliabilityCancelaPetingi(0, 4, 4);

    EXPECT_NEAR(result1.reliability, 0.9997948026, 1e-10);
    EXPECT_EQ(result1.recursions, result2.recursions) << "Recursion count must be deterministic";
    EXPECT_EQ(result1.recursions, 206);  // Optimized with corrected ISPT
}

TEST_F(ReliabilityGraphTest, EmptyGraph) {
    ReliabilityGraph empty({0}, {}, {}, {});
    
    EXPECT_EQ(empty.numVertices(), 0);
    EXPECT_EQ(empty.numEdges(), 0);
    EXPECT_TRUE(empty.target_vertices_.empty());
}

TEST_F(ReliabilityGraphTest, SingleVertexGraph) {
    std::vector<ReliabilityGraph::VertexId> kao = {0, 0};
    std::vector<ReliabilityGraph::VertexId> fo = {};
    std::vector<ReliabilityGraph::Probability> p = {};
    std::vector<int> targets = {1};
    
    ReliabilityGraph graph(kao, fo, p, targets);
    
    EXPECT_EQ(graph.numVertices(), 1);
    EXPECT_EQ(graph.numEdges(), 0);
    
    // Reliability from vertex to itself should be 1.0
    auto result = graph.calculateReliabilityBetweenVertices(0, 0, 1);
    EXPECT_EQ(result.reliability, 1.0);
}
