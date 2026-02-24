/**
 * @file Logger.cpp
 * @brief Implementation of the Logger class
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/Logger.h"

#ifdef GRAPH_RELIABILITY_ENABLE_LOGGING

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <vector>
#include <memory>

namespace graph_reliability {

std::shared_ptr<spdlog::logger> Logger::logger_ = nullptr;
bool Logger::initialized_ = false;

void Logger::initialize(const std::string& log_file_path,
                       size_t max_file_size,
                       size_t max_files) {
    if (initialized_) {
        return;
    }

    try {
        // Create sinks
        std::vector<spdlog::sink_ptr> sinks;

        // Console sink with colors
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::info);
        console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
        sinks.push_back(console_sink);

        // File sink with rotation
        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            log_file_path, max_file_size, max_files);
        file_sink->set_level(spdlog::level::debug);
        file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] [%t] %v");
        sinks.push_back(file_sink);

        // Create logger with both sinks
        logger_ = std::make_shared<spdlog::logger>("graph_reliability", sinks.begin(), sinks.end());
        logger_->set_level(spdlog::level::debug);
        logger_->flush_on(spdlog::level::warn);

        // Register as default logger
        spdlog::register_logger(logger_);
        spdlog::set_default_logger(logger_);

        initialized_ = true;
    } catch (const spdlog::spdlog_ex& ex) {
        // If logging initialization fails, we can't log it, so just set to null
        logger_ = nullptr;
        initialized_ = false;
    }
}

void Logger::shutdown() {
    if (logger_) {
        logger_->flush();
        spdlog::shutdown();
        logger_ = nullptr;
        initialized_ = false;
    }
}

void Logger::setLevel(const std::string& level) {
    if (!logger_) {
        return;
    }

    if (level == "trace") {
        logger_->set_level(spdlog::level::trace);
    } else if (level == "debug") {
        logger_->set_level(spdlog::level::debug);
    } else if (level == "info") {
        logger_->set_level(spdlog::level::info);
    } else if (level == "warn") {
        logger_->set_level(spdlog::level::warn);
    } else if (level == "error") {
        logger_->set_level(spdlog::level::err);
    } else if (level == "critical") {
        logger_->set_level(spdlog::level::critical);
    } else if (level == "off") {
        logger_->set_level(spdlog::level::off);
    }
}

std::shared_ptr<spdlog::logger> Logger::getLogger() {
    if (!initialized_ && !logger_) {
        // Auto-initialize with defaults if not already initialized
        initialize();
    }
    return logger_;
}

} // namespace graph_reliability

#endif // GRAPH_RELIABILITY_ENABLE_LOGGING
