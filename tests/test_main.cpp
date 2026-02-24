/**
 * @file test_main.cpp
 * @brief Main entry point for unit tests
 * @author Graduate Work Project
 * @date 2026
 */

#include <gtest/gtest.h>
#include "graph_reliability/Logger.h"

int main(int argc, char **argv) {
    // Initialize logger for tests (if logging is enabled)
    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    graph_reliability::Logger::initialize("test.log");
    graph_reliability::Logger::setLevel("warn"); // Reduce test output noise
    #endif

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    graph_reliability::Logger::shutdown();
    #endif

    return result;
}
