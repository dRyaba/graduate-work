/**
 * @file test_Graph.cpp
 * @brief Unit tests for Graph class
 * @author Graduate Work Project
 * @date 2026
 */

#include <gtest/gtest.h>
#include "graph_reliability/Graph.h"
#include "graph_reliability/Exceptions.h"

using namespace graph_reliability;

class GraphTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple graph: triangle (3 vertices, 3 edges)
        // Vertices: 0, 1, 2
        // Edges: 0-1, 1-2, 2-0
        kao_ = {0, 2, 4, 6};
        fo_ = {1, 2, 0, 2, 0, 1};
        p_array_ = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    }

    std::vector<Graph::VertexId> kao_;
    std::vector<Graph::VertexId> fo_;
    std::vector<Graph::Probability> p_array_;
};

TEST_F(GraphTest, ConstructorAndValidation) {
    Graph graph(kao_, fo_, p_array_);
    
    EXPECT_EQ(graph.numVertices(), 3);
    EXPECT_EQ(graph.numEdges(), 6); // Undirected: 3 edges * 2 directions
}

TEST_F(GraphTest, InvalidGraphThrowsException) {
    // Invalid: KAO doesn't start with 0
    std::vector<Graph::VertexId> invalid_kao = {1, 2, 4, 6};
    EXPECT_THROW(Graph(invalid_kao, fo_, p_array_), InvalidGraphException);
    
    // Invalid: KAO last element doesn't match FO size
    std::vector<Graph::VertexId> invalid_kao2 = {0, 2, 4, 10};
    EXPECT_THROW(Graph(invalid_kao2, fo_, p_array_), InvalidGraphException);
    
    // Invalid: FO and p_array size mismatch
    std::vector<Graph::Probability> invalid_p = {0.9, 0.9};
    EXPECT_THROW(Graph(kao_, fo_, invalid_p), InvalidGraphException);
}

TEST_F(GraphTest, NumVerticesAndEdges) {
    Graph graph(kao_, fo_, p_array_);
    
    EXPECT_EQ(graph.numVertices(), 3);
    EXPECT_EQ(graph.numEdges(), 6);
    
    // Empty graph
    Graph empty({0}, {}, {});
    EXPECT_EQ(empty.numVertices(), 0);
    EXPECT_EQ(empty.numEdges(), 0);
}

TEST_F(GraphTest, HasEdge) {
    Graph graph(kao_, fo_, p_array_);
    
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(1, 0));
    EXPECT_TRUE(graph.hasEdge(1, 2));
    EXPECT_TRUE(graph.hasEdge(2, 0));
    EXPECT_FALSE(graph.hasEdge(0, 3)); // Out of range
    EXPECT_FALSE(graph.hasEdge(3, 0)); // Out of range
}

TEST_F(GraphTest, FindEdge) {
    Graph graph(kao_, fo_, p_array_);
    
    auto edge_id = graph.findEdge(0, 1);
    EXPECT_GE(edge_id, 0);
    EXPECT_LT(edge_id, graph.numEdges());
    
    auto invalid_edge = graph.findEdge(0, 3);
    EXPECT_EQ(invalid_edge, -1);
}

TEST_F(GraphTest, RemoveEdge) {
    Graph graph(kao_, fo_, p_array_);
    
    Graph modified = graph.removeEdge(0, 1);
    
    EXPECT_FALSE(modified.hasEdge(0, 1));
    EXPECT_FALSE(modified.hasEdge(1, 0));
    EXPECT_EQ(modified.numEdges(), graph.numEdges() - 2); // Undirected edge removed
}

TEST_F(GraphTest, CalculateDistance) {
    Graph graph(kao_, fo_, p_array_);
    
    // Distance from vertex to itself
    EXPECT_EQ(graph.calculateDistance(0, 0, 10), 0);
    
    // Distance between connected vertices
    EXPECT_EQ(graph.calculateDistance(0, 1, 10), 1);
    
    // Distance in triangle
    EXPECT_EQ(graph.calculateDistance(0, 2, 10), 1);
    
    // Max distance exceeded
    int max_dist = 0;
    EXPECT_GT(graph.calculateDistance(0, 1, max_dist), max_dist);
}

TEST_F(GraphTest, SwapVertices) {
    Graph graph(kao_, fo_, p_array_);
    
    Graph swapped = graph.swapVertices(0, 1);
    
    // After swapping 0 and 1, edge (0,1) should still exist but vertices are swapped
    EXPECT_TRUE(swapped.hasEdge(1, 0)); // Was (0,1), now (1,0)
    EXPECT_EQ(swapped.numVertices(), graph.numVertices());
    EXPECT_EQ(swapped.numEdges(), graph.numEdges());
}

TEST_F(GraphTest, EmptyGraph) {
    Graph empty({0}, {}, {});
    
    EXPECT_EQ(empty.numVertices(), 0);
    EXPECT_EQ(empty.numEdges(), 0);
    EXPECT_FALSE(empty.hasEdge(0, 1));
}
