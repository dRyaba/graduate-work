/**
 * @file CliHandlers.h
 * @brief Command dispatch handlers for the graph_reliability CLI.
 *
 * Each handler implements one top-level subcommand (args[0] is the flag name,
 * e.g. "--run"). The main() entry point parses the global "--verbose" flag
 * and routes the remaining arguments to the appropriate handler.
 */

#pragma once

#include <string>
#include <vector>

namespace graph_reliability {

/** Print full CLI usage to stdout. */
void printUsage();

/** @brief Handle `--run` single-graph test. */
int handleRun(const std::vector<std::string>& args);

/** @brief Handle `--convert edge2kao | kao2edge` format conversions. */
int handleConversion(const std::vector<std::string>& args);

/** @brief Handle `--test <method_id>` comprehensive test suite. */
int handleTesting(const std::vector<std::string>& args);

/** @brief Handle `--cross-check` all-methods consistency run. */
int handleCrossCheck(const std::vector<std::string>& args);

/** @brief Handle `--visualize <graph> <out>` layout + SVG/DOT export. */
int handleVisualize(const std::vector<std::string>& args);

} // namespace graph_reliability
