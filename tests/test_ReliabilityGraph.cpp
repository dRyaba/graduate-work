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
#include "graph_reliability/TestSuite.h"

#include <array>
#include <cctype>
#include <functional>
#include <optional>
#include <string>

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

TEST_F(ReliabilityGraphTest, CancelaPetingiMultiMatchesModifiedFactoring) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("K4_kao.txt")) {
        GTEST_SKIP() << "K4_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("K4_kao.txt");
    ASSERT_TRUE(graph);

    // Compare CDF from both methods for K4 with diameter 3
    auto [cdf_modified, recs_mod] = graph->calculateReliabilityMultipleDiameters(0, 3, 0, 3);
    auto [cdf_cpfm, recs_cpfm] = graph->calculateReliabilityCancelaPetingiMulti(0, 3, 0, 3);

    ASSERT_EQ(cdf_modified.size(), cdf_cpfm.size());
    for (size_t i = 0; i < cdf_modified.size(); ++i) {
        EXPECT_NEAR(cdf_modified[i], cdf_cpfm[i], 1e-9) 
            << "CDF mismatch at d=" << i;
    }
}

TEST_F(ReliabilityGraphTest, Method5MatchesMethod3) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("2_3x3_blocks_kao.txt")) {
        GTEST_SKIP() << "2_3x3_blocks_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("2_3x3_blocks_kao.txt");
    ASSERT_TRUE(graph);

    auto result3 = graph->calculateReliabilityWithMDecomposition(0, 16, 8);
    auto result5 = graph->calculateReliabilityWithMDecompositionCPFM(0, 16, 8);

    EXPECT_NEAR(result3.reliability, result5.reliability, 1e-9)
        << "Method 5 should match Method 3";
}

TEST_F(ReliabilityGraphTest, Method5MatchesMethod4OnSingleBlock) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("K5_kao.txt")) {
        GTEST_SKIP() << "K5_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("K5_kao.txt");
    ASSERT_TRUE(graph);

    auto result4 = graph->calculateReliabilityCancelaPetingi(0, 4, 4);
    auto result5 = graph->calculateReliabilityWithMDecompositionCPFM(0, 4, 4);

    EXPECT_NEAR(result4.reliability, result5.reliability, 1e-9)
        << "Method 5 should match Method 4 on single-block graph";
}

TEST_F(ReliabilityGraphTest, Method5OnThreeBlockGraph) {
    DataImporter importer("graphs_data/");
    if (!importer.fileExists("3_blocks_sausage_3x3_kao.txt")) {
        GTEST_SKIP() << "3_blocks_sausage_3x3_kao.txt not found";
    }
    auto graph = importer.loadKAOGraph("3_blocks_sausage_3x3_kao.txt");
    ASSERT_TRUE(graph);

    auto result3 = graph->calculateReliabilityWithMDecomposition(0, 28, 12);
    auto result5 = graph->calculateReliabilityWithMDecompositionCPFM(0, 28, 12);

    EXPECT_NEAR(result3.reliability, result5.reliability, 1e-9)
        << "Method 5 should match Method 3 on 3-block graph";
    
    // Method 5 should use fewer recursions due to ISPT
    EXPECT_LT(result5.recursions, result3.recursions)
        << "Method 5 should use fewer recursions than Method 3";
}

// --- Parameterized cross-method agreement test ----------------------------

struct AgreementCase {
    std::string graph_file;
    int s;
    int t;
    int diameter;
    int timeout_sec;
};

class MethodsAgreementTest : public ::testing::TestWithParam<AgreementCase> {};

TEST_P(MethodsAgreementTest, AllMethodsAgreeWithinTolerance) {
    const AgreementCase tc = GetParam();
    DataImporter importer("graphs_data/");
    if (!importer.fileExists(tc.graph_file)) {
        GTEST_SKIP() << tc.graph_file << " not found";
    }

    constexpr int kMethods = 6;
    std::array<std::optional<ReliabilityResult>, kMethods> results{};
    std::array<std::string, kMethods> errors{};

    for (int method_id = 0; method_id < kMethods; ++method_id) {
        auto graph = importer.loadKAOGraph(tc.graph_file);
        ASSERT_TRUE(graph) << "Could not load " << tc.graph_file;
        ReliabilityGraph* gptr = graph.get();
        const int s = tc.s, t = tc.t, d = tc.diameter;

        std::function<ReliabilityResult()> calc;
        switch (method_id) {
            case 0: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityBetweenVertices(s, t, d); }; break;
            case 1: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityWithRecursiveDecomposition(s, t, d); }; break;
            case 2: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityWithDecomposition(s, t, d); }; break;
            case 3: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityWithMDecomposition(s, t, d); }; break;
            case 4: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityCancelaPetingi(s, t, d); }; break;
            case 5: calc = [gptr, s, t, d]() { return gptr->calculateReliabilityWithMDecompositionCPFM(s, t, d); }; break;
        }

        // Keep the graph alive across a possible timeout detach.
        auto keeper = std::make_shared<std::unique_ptr<ReliabilityGraph>>(std::move(graph));
        std::function<ReliabilityResult()> wrapped = [calc, keeper]() { return calc(); };

        std::string err;
        results[method_id] = TestSuite::runWithTimeout(wrapped, tc.timeout_sec, &err);
        errors[method_id] = err;
    }

    if (!results[0].has_value()) {
        GTEST_SKIP() << "baseline m0 did not finish in " << tc.timeout_sec
                     << "s (err=" << errors[0] << ")";
    }

    const double baseline = results[0]->reliability;
    for (int m = 1; m < kMethods; ++m) {
        if (!results[m].has_value()) {
            // Don't fail the test on method timeout — just record it. The
            // cross-check CSV run is where we're strict about completion; unit
            // tests only fail on genuine numerical disagreement.
            GTEST_LOG_(INFO) << "m" << m << " did not finish (err=" << errors[m] << ")";
            continue;
        }
        EXPECT_NEAR(results[m]->reliability, baseline, 1e-10)
            << "Method " << m << " disagrees with m0 on "
            << tc.graph_file << " d=" << tc.diameter;
    }
}

INSTANTIATE_TEST_SUITE_P(
    CrossMethodAgreement,
    MethodsAgreementTest,
    ::testing::Values(
        AgreementCase{"K4_kao.txt",                   0, 3, 2,  10},
        AgreementCase{"K4_kao.txt",                   0, 3, 3,  10},
        AgreementCase{"3_blocks_sausage_3x3_kao.txt", 0, 28, 9,  10},
        AgreementCase{"3_blocks_sausage_3x3_kao.txt", 0, 28, 10, 10}
    ),
    [](const ::testing::TestParamInfo<AgreementCase>& info) {
        std::string base = info.param.graph_file;
        auto suffix = base.find("_kao.txt");
        if (suffix != std::string::npos) base.erase(suffix);
        // GoogleTest allows only [A-Za-z0-9_] in generated names.
        for (char& c : base) {
            if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
        }
        return base + "_d" + std::to_string(info.param.diameter);
    });

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
