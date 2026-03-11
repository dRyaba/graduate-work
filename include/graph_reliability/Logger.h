/**
 * @file Logger.h
 * @brief Logging system wrapper with conditional compilation for zero overhead in Release
 * @author Graduate Work Project
 * @date 2026
 * 
 * This header provides a logging interface that completely removes logging code
 * in Release builds through conditional compilation. In Debug builds, it uses spdlog
 * for efficient logging with file rotation and console output.
 */

#pragma once

#ifdef GRAPH_RELIABILITY_ENABLE_LOGGING

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <memory>
#include <string>

namespace graph_reliability {

/**
 * @brief Logger class wrapper around spdlog
 * 
 * Provides logging functionality that is completely disabled in Release builds.
 * In Debug builds, supports both console and file logging with rotation.
 */
class Logger {
public:
    /**
     * @brief Initialize logger with console and file sinks
     * @param log_file_path Path to log file (default: "graph_reliability.log")
     * @param max_file_size Maximum log file size in bytes (default: 5MB)
     * @param max_files Maximum number of rotated files (default: 3)
     */
    static void initialize(const std::string& log_file_path = "graph_reliability.log",
                          size_t max_file_size = 5 * 1024 * 1024,
                          size_t max_files = 3);

    /**
     * @brief Shutdown logger
     */
    static void shutdown();

    /**
     * @brief Set logging level
     * @param level Logging level (trace, debug, info, warn, error, critical, off)
     */
    static void setLevel(const std::string& level);

    /**
     * @brief Get the underlying spdlog logger instance
     * @return Shared pointer to logger
     */
    static std::shared_ptr<spdlog::logger> getLogger();

private:
    static std::shared_ptr<spdlog::logger> logger_;
    static bool initialized_;
};

} // namespace graph_reliability

// Macro definitions for logging - these will be empty in Release builds
#define LOG_TRACE(...) graph_reliability::Logger::getLogger()->trace(__VA_ARGS__)
#define LOG_DEBUG(...) graph_reliability::Logger::getLogger()->debug(__VA_ARGS__)
#define LOG_INFO(...) graph_reliability::Logger::getLogger()->info(__VA_ARGS__)
#define LOG_WARN(...) graph_reliability::Logger::getLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...) graph_reliability::Logger::getLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...) graph_reliability::Logger::getLogger()->critical(__VA_ARGS__)

#else // GRAPH_RELIABILITY_ENABLE_LOGGING not defined

#include <cstddef>
#include <string>

// Empty macros - compiler will completely remove these calls in Release builds
#define LOG_TRACE(...) ((void)0)
#define LOG_DEBUG(...) ((void)0)
#define LOG_INFO(...) ((void)0)
#define LOG_WARN(...) ((void)0)
#define LOG_ERROR(...) ((void)0)
#define LOG_CRITICAL(...) ((void)0)

// Dummy Logger class for Release builds (for code that might reference it)
namespace graph_reliability {
class Logger {
public:
    static void initialize(const std::string& = "", size_t = 0, size_t = 0) {}
    static void shutdown() {}
    static void setLevel(const std::string&) {}
};
} // namespace graph_reliability

#endif // GRAPH_RELIABILITY_ENABLE_LOGGING
