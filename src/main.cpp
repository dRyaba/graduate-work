/**
 * @file main.cpp
 * @brief Entry point for the Graph Reliability Analysis CLI.
 *
 * Parses global flags (--verbose) and dispatches to the handler for the
 * chosen subcommand. All subcommand logic lives in CliHandlers.cpp.
 */

#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include <cstdio>

#include "graph_reliability.h"
#include "graph_reliability/Logger.h"
#include "graph_reliability/CliHandlers.h"

using namespace graph_reliability;

int main(int argc, char* argv[]) {
    // Flush stdout/stderr on every write so progress survives a redirect
    // to file and is not lost on TaskStop/Ctrl+C. MinGW's C stdio defaults
    // to full-buffering when stdout is not a terminal; on top of that
    // std::cout has its own buffer. We address both: setvbuf covers the
    // C FILE* path (printf/LOG) and std::unitbuf makes std::cout flush
    // after every << operation.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    std::setvbuf(stderr, nullptr, _IOLBF, 0);
    std::cout << std::unitbuf;
    std::cerr << std::unitbuf;

    #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
    bool verbose = false;
    std::vector<std::string> args(argv + 1, argv + argc);
    for (const auto& arg : args) {
        if (arg == "--verbose" || arg == "-v") {
            verbose = true;
            break;
        }
    }

    Logger::initialize();
    Logger::setLevel(verbose ? "debug" : "info");
    LOG_INFO("Graph Reliability Analysis Tool v{} started", getVersion());
    #else
    std::vector<std::string> args(argv + 1, argv + argc);
    #endif

    if (args.empty() || args[0] == "--help") {
        printUsage();
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        Logger::shutdown();
        #endif
        return 0;
    }

    try {
        int result = 0;
        if (args[0] == "--convert") {
            result = handleConversion(args);
        } else if (args[0] == "--test") {
            result = handleTesting(args);
        } else if (args[0] == "--cross-check") {
            result = handleCrossCheck(args);
        } else if (args[0] == "--run") {
            result = handleRun(args);
        } else if (args[0] == "--visualize") {
            result = handleVisualize(args);
        } else {
            LOG_ERROR("Unknown command: {}", args[0]);
            std::cerr << "Error: Unknown command: " << args[0] << "\n";
            std::cerr << "Use --help for usage information.\n";
            result = 1;
        }

        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        LOG_INFO("Application exiting with code {}", result);
        Logger::shutdown();
        #endif

        return result;
    } catch (const std::exception& e) {
        LOG_CRITICAL("Fatal error: {}", e.what());
        std::cerr << "Fatal error: " << e.what() << "\n";
        #ifdef GRAPH_RELIABILITY_ENABLE_LOGGING
        Logger::shutdown();
        #endif
        return 2;
    }
}
