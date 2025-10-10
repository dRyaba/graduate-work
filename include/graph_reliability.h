/**
 * @file graph_reliability.h
 * @brief Main header file for Graph Reliability Analysis Library
 * @author Graduate Work Project
 * @date 2024
 * 
 * This library provides comprehensive tools for network reliability analysis,
 * including graph representation, reliability calculations, and testing frameworks.
 */

#pragma once

// Core graph structures
#include "graph_reliability/Graph.h"
#include "graph_reliability/ReliabilityGraph.h"

// Data management
#include "graph_reliability/DataImporter.h"

// Testing framework
#include "graph_reliability/TestSuite.h"

// Utilities
#include "graph_reliability/GraphOperations.h"

/**
 * @namespace graph_reliability
 * @brief Main namespace for the Graph Reliability Analysis Library
 * 
 * This namespace contains all classes and functions related to graph reliability
 * analysis, providing a clean and organized API for network reliability calculations.
 */
namespace graph_reliability {

/**
 * @brief Library version information
 */
struct Version {
    static constexpr int major = 1;
    static constexpr int minor = 0;
    static constexpr int patch = 0;
    static constexpr const char* string = "1.0.0";
};

/**
 * @brief Get library version string
 * @return Version string
 */
inline const char* getVersion() noexcept {
    return Version::string;
}

/**
 * @brief Get library version components
 * @return Version struct
 */
inline Version getVersionInfo() noexcept {
    return Version{};
}

} // namespace graph_reliability
