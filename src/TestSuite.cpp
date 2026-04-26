/**
 * @file TestSuite.cpp
 * @brief Implementation of the TestSuite class
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/TestSuite.h"
#include "graph_reliability/GraphOperations.h"
#include "graph_reliability/Logger.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <future>
#include <thread>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <array>

namespace graph_reliability {

TestSuite::TestSuite(DataImporter& data_importer)
    : data_importer_(data_importer), number_of_repetitions_(1) {
}

void TestSuite::setTestConfigurations(const std::map<std::string, std::vector<int>>& configs) {
    test_configurations_ = configs;
}

std::vector<std::string> TestSuite::getAvailableMethods() const {
    return {
        "Standard Factoring (Level 0)",
        "Recursive Decomposition (Level 1)",
        "Simple Factoring (Level 2)",
        "M-Decomposition (Level 3)",
        "Cancela-Petingi (Level 4)",
        "M-Decomp + CPFM (Level 5)"
    };
}

TestRunResult TestSuite::runSingleTest(const TestConfiguration& config, int method_id) {
    std::vector<double> times;
    std::vector<double> reliabilities;
    std::vector<long long> recursions;

    LOG_INFO("Running single test: file={}, diameter={}, reps={}, method={}",
             config.graph_filename, config.upper_bound_diameter, 
             config.number_of_repetitions, method_id);
    std::cout << "Testing " << config.graph_filename 
              << " (D=" << config.upper_bound_diameter << ")"
              << " Reps=" << config.number_of_repetitions << "..." << std::endl;

    for (int i = 0; i < config.number_of_repetitions; ++i) {
        try {
            LOG_DEBUG("Test repetition {}/{}", i + 1, config.number_of_repetitions);
            // Load fresh graph for each run to ensure clean state
            auto graph = data_importer_.loadKAOGraph(config.graph_filename);
            
            if (!graph) {
                LOG_ERROR("Failed to load graph: {}", config.graph_filename);
                std::cerr << "Failed to load graph: " << config.graph_filename << std::endl;
                continue;
            }

            int s = config.source_vertex;
            int t = config.target_vertex;

            ReliabilityResult result;
            
            switch (method_id) {
                case 0: // Level 0: Standard Factoring (Baseline - SLOWEST)
                    LOG_DEBUG("Using Standard Factoring method");
                    result = graph->calculateReliabilityBetweenVertices(s, t, config.upper_bound_diameter);
                    break;
                case 1: // Level 1: Recursive Decomposition (Nested Recursion - INEFFICIENT)
                    LOG_DEBUG("Using Recursive Decomposition method");
                    result = graph->calculateReliabilityWithRecursiveDecomposition(s, t, config.upper_bound_diameter);
                    break;
                case 2: // Level 2: Simple Factoring (Convolution - FAST)
                    LOG_DEBUG("Using Simple Factoring method");
                    result = graph->calculateReliabilityWithDecomposition(s, t, config.upper_bound_diameter);
                    break;
                case 3: // Level 3: M-Decomposition (FASTEST)
                    LOG_DEBUG("Using M-Decomposition method");
                    result = graph->calculateReliabilityWithMDecomposition(s, t, config.upper_bound_diameter);
                    break;
                case 4: // Level 4: Cancela-Petingi path-based factoring
                    LOG_DEBUG("Using Cancela-Petingi method");
                    result = graph->calculateReliabilityCancelaPetingi(s, t, config.upper_bound_diameter);
                    break;
                case 5: // Level 5: M-Decomposition + CPFM hybrid
                    LOG_DEBUG("Using M-Decomp + CPFM method");
                    result = graph->calculateReliabilityWithMDecompositionCPFM(s, t, config.upper_bound_diameter);
                    break;
                default:
                    LOG_ERROR("Unknown method ID: {}", method_id);
                    std::cerr << "Unknown method ID: " << method_id << std::endl;
                    return TestRunResult();
            }

            LOG_DEBUG("Repetition {} result: reliability={}, time={}s, recursions={}",
                      i + 1, result.reliability, result.execution_time_sec, result.recursions);

            times.push_back(result.execution_time_sec);
            reliabilities.push_back(result.reliability);
            recursions.push_back(result.recursions);

        } catch (const std::exception& e) {
            LOG_ERROR("Error in test run {}: {}", i + 1, e.what());
            std::cerr << "Error in test run " << i + 1 << ": " << e.what() << std::endl;
        }
    }

    if (times.empty()) return TestRunResult();

    double avg_time = 0;
    double avg_rel = 0;
    long long avg_rec = 0;

    for (double t : times) avg_time += t;
    for (double r : reliabilities) avg_rel += r;
    for (long long rec : recursions) avg_rec += rec;

    TestRunResult result(
        avg_time / times.size(),
        avg_rel / reliabilities.size(),
        avg_rec / recursions.size()
    );
    
    LOG_INFO("Test completed: avg_reliability={}, avg_time={}s, avg_recursions={}",
             result.average_reliability, result.average_time_seconds, result.average_recursions);
    
    return result;
}

bool TestSuite::runComprehensiveTests(int method_id, const std::string& output_filename) {
    LOG_INFO("Starting comprehensive tests: method_id={}, output_file={}", method_id, output_filename);
    std::vector<std::pair<TestConfiguration, TestRunResult>> results;
    bool all_success = true;
    std::vector<std::string> methods = getAvailableMethods();
    std::string method_name = (method_id >= 0 && method_id < methods.size()) ? methods[method_id] : "Unknown";

    LOG_INFO("Test configurations: {} graphs", test_configurations_.size());
    for (const auto& [filename, diameters] : test_configurations_) {
        try {
            LOG_DEBUG("Processing graph: {}, diameters: {}", filename, diameters.size());
            // Find s/t for this graph
            auto temp_graph = data_importer_.loadKAOGraph(filename);
            if (!temp_graph) {
                LOG_ERROR("Skipping {}: could not load", filename);
                std::cerr << "Skipping " << filename << ": could not load" << std::endl;
                all_success = false;
                continue;
            }

            auto [s, t] = findSourceAndTargetVertices(*temp_graph);
            if (s == -1 || t == -1) {
                LOG_ERROR("Skipping {}: could not find source/target vertices", filename);
                std::cerr << "Skipping " << filename << ": could not find source/target vertices" << std::endl;
                all_success = false;
                continue;
            }

            LOG_DEBUG("Found source={}, target={} for graph {}", s, t, filename);

            for (int d : diameters) {
                TestConfiguration config(filename, s, t, d, number_of_repetitions_);
                TestRunResult result = runSingleTest(config, method_id);
                results.push_back({config, result});
            }

        } catch (const std::exception& e) {
            LOG_ERROR("Error processing {}: {}", filename, e.what());
            std::cerr << "Error processing " << filename << ": " << e.what() << std::endl;
            all_success = false;
        }
    }

    LOG_INFO("Writing {} test results to {}", results.size(), output_filename);
    writeResultsToCSV(results, method_name, output_filename);
    LOG_INFO("Comprehensive tests completed: success={}", all_success);
    return all_success;
}

std::pair<int, int> TestSuite::findSourceAndTargetVertices(const ReliabilityGraph& graph) const {
    // Logic from legacy main.cpp:
    // Find vertices where Targets[i] == 1. First is s, second is t.
    // Note: ReliabilityGraph uses 0-based target_vertices_ vector.
    // But legacy Targets was 1-based vector (size N+1).
    // My DataImporter converts Targets to 0-based vector of size N?
    // Or size N+1? 
    // Let's check DataImporter.cpp: 
    // "if (!targets.empty()) targets.erase(targets.begin());" -> Size N.
    // So indices are 0..N-1.

    int s = -1, t = -1;
    const auto& targets = graph.target_vertices_;
    
    for (size_t i = 0; i < targets.size(); ++i) {
        if (targets[i] == 1) {
            if (s == -1) s = static_cast<int>(i);
            else if (t == -1) t = static_cast<int>(i);
        }
    }
    
    return {s, t};
}

void TestSuite::writeResultsToCSV(const std::vector<std::pair<TestConfiguration, TestRunResult>>& results,
                                 const std::string& method_name,
                                 const std::string& filename) const {
    LOG_DEBUG("Writing CSV results to: {}", filename);
    std::ofstream file(filename);
    if (!file) {
        LOG_ERROR("Could not open output file: {}", filename);
        std::cerr << "Could not open output file: " << filename << std::endl;
        return;
    }

    file << "Graph,S,T,D,Reps,Method,AvgTime,AvgRel,AvgRecs\n";

    for (const auto& [config, result] : results) {
        file << config.graph_filename << ","
             << config.source_vertex << ","
             << config.target_vertex << ","
             << config.upper_bound_diameter << ","
             << config.number_of_repetitions << ","
             << method_name << ","
             << std::fixed << std::setprecision(6) << result.average_time_seconds << ","
             << std::setprecision(15) << result.average_reliability << ","
             << result.average_recursions << "\n";
    }

    LOG_INFO("Successfully wrote {} results to {}", results.size(), filename);
    std::cout << "Results written to " << filename << std::endl;
}

// ---- Cross-check implementation -------------------------------------------

std::optional<ReliabilityResult> TestSuite::runWithTimeout(
        std::function<ReliabilityResult()> calc,
        int timeout_sec,
        std::string* out_error) {
    // std::async(launch::async) guarantees a new thread. We keep the future
    // on the heap so that on timeout we can detach (leak) and return without
    // the future's destructor blocking on join. The worker thread then runs
    // until the process exits — acceptable for a one-shot cross-check.
    struct Shared {
        std::promise<ReliabilityResult> promise;
        std::atomic<bool> claimed{false};   // who sets the promise first wins
    };
    auto shared = std::make_shared<Shared>();
    auto fut = shared->promise.get_future();

    std::thread worker([shared, calc = std::move(calc)]() {
        try {
            ReliabilityResult r = calc();
            if (!shared->claimed.exchange(true)) {
                shared->promise.set_value(r);
            }
        } catch (...) {
            if (!shared->claimed.exchange(true)) {
                try { shared->promise.set_exception(std::current_exception()); }
                catch (...) { /* nothing to do */ }
            }
        }
    });

    auto status = fut.wait_for(std::chrono::seconds(std::max(1, timeout_sec)));
    if (status == std::future_status::ready) {
        worker.join();
        try {
            return fut.get();
        } catch (const std::exception& e) {
            if (out_error) *out_error = e.what();
            return std::nullopt;
        } catch (...) {
            if (out_error) *out_error = "unknown exception";
            return std::nullopt;
        }
    }

    // Timeout: detach the worker and move on. Block any future set_value().
    shared->claimed.store(true);
    worker.detach();
    return std::nullopt;
}

std::optional<ReliabilityCdfResult> TestSuite::runWithTimeoutCdf(
        std::function<ReliabilityCdfResult()> calc,
        int timeout_sec,
        std::string* out_error) {
    struct Shared {
        std::promise<ReliabilityCdfResult> promise;
        std::atomic<bool> claimed{false};
    };
    auto shared = std::make_shared<Shared>();
    auto fut = shared->promise.get_future();

    std::thread worker([shared, calc = std::move(calc)]() {
        try {
            ReliabilityCdfResult r = calc();
            if (!shared->claimed.exchange(true)) {
                shared->promise.set_value(std::move(r));
            }
        } catch (...) {
            if (!shared->claimed.exchange(true)) {
                try { shared->promise.set_exception(std::current_exception()); }
                catch (...) {}
            }
        }
    });

    auto status = fut.wait_for(std::chrono::seconds(std::max(1, timeout_sec)));
    if (status == std::future_status::ready) {
        worker.join();
        try {
            return fut.get();
        } catch (const std::exception& e) {
            if (out_error) *out_error = e.what();
            return std::nullopt;
        } catch (...) {
            if (out_error) *out_error = "unknown exception";
            return std::nullopt;
        }
    }

    shared->claimed.store(true);
    worker.detach();
    return std::nullopt;
}

std::vector<int> TestSuite::chooseDiameters(const ReliabilityGraph& graph,
                                            int s, int t,
                                            int count, int step) {
    if (count <= 0 || step <= 0) return {};
    const int n = static_cast<int>(graph.target_vertices_.size());
    if (n <= 1 || s < 0 || t < 0 || s >= n || t >= n) return {};

    const int max_dist = n - 1;
    int dist = graph.calculateDistance(s, t, max_dist);
    if (dist <= 0 || dist > max_dist) return {};  // unreachable or self-loop

    std::vector<int> out;
    out.reserve(count);
    for (int i = 0; i < count; ++i) {
        int d = dist + i * step;
        if (d > max_dist) break;
        out.push_back(d);
    }
    return out;
}

static const char* statusToString(CrossCheckStatus s) {
    switch (s) {
        case CrossCheckStatus::OK:      return "OK";
        case CrossCheckStatus::TIMEOUT: return "TIMEOUT";
        case CrossCheckStatus::ERROR:   return "ERROR";
    }
    return "UNKNOWN";
}

void TestSuite::writeCrossCheckCSV(const std::vector<CrossCheckRow>& rows,
                                   const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Could not open output file: " << filename << std::endl;
        return;
    }
    file << "Graph,S,T,D,MethodId,Method,Status,Reliability,TimeSec,Recursions,DiffFromBaseline,Error\n";
    // Group by (graph, d) to compute DiffFromBaseline.
    // Baseline = first successful method in order: 0, 3, 5, then any OK.
    auto key = [](const CrossCheckRow& r) {
        return r.graph + "|" + std::to_string(r.diameter);
    };
    std::map<std::string, double> baseline;
    for (int preferred : {0, 3, 5}) {
        for (const auto& r : rows) {
            if (r.method_id != preferred) continue;
            if (r.status != CrossCheckStatus::OK) continue;
            baseline.emplace(key(r), r.reliability);
        }
    }
    for (const auto& r : rows) {
        if (r.status != CrossCheckStatus::OK) continue;
        baseline.emplace(key(r), r.reliability);
    }

    for (const auto& r : rows) {
        std::string diff_str;
        if (r.status == CrossCheckStatus::OK) {
            auto it = baseline.find(key(r));
            if (it != baseline.end()) {
                std::ostringstream oss;
                oss << std::scientific << std::setprecision(3)
                    << std::abs(r.reliability - it->second);
                diff_str = oss.str();
            }
        }
        std::string rel_str;
        if (r.status == CrossCheckStatus::OK) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(15) << r.reliability;
            rel_str = oss.str();
        }
        std::string time_str;
        if (r.status != CrossCheckStatus::TIMEOUT) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(6) << r.time_seconds;
            time_str = oss.str();
        }
        std::string err = r.error_message;
        // CSV-escape: replace comma/quote with space for simplicity.
        std::replace(err.begin(), err.end(), ',', ' ');
        std::replace(err.begin(), err.end(), '"', ' ');
        file << r.graph << ","
             << r.s << ","
             << r.t << ","
             << r.diameter << ","
             << r.method_id << ","
             << r.method_name << ","
             << statusToString(r.status) << ","
             << rel_str << ","
             << time_str << ","
             << r.recursions << ","
             << diff_str << ","
             << err << "\n";
    }
    std::cout << "Cross-check results written to " << filename << std::endl;
}

namespace {

// Per-method timeouts in seconds. Picked to match observed difficulty:
// m0 is the reference but gets enough budget to clear sausage-3 d=10; m1 is
// known-slow and we cap it tight; m3/m5 are the workhorses and get the most
// budget; m2/m4 sit in the middle.
constexpr std::array<int, 6> kMethodTimeoutsSec = {
    60,  // m0: Standard Factoring
    30,  // m1: Recursive Decomposition (known slow)
    120, // m2: Simple Factoring
    300, // m3: M-Decomposition (multi-d in cross-check)
    300, // m4: Cancela-Petingi (multi-d enumerates all paths up to d_max)
    300, // m5: M-Decomp + CPFM (multi-d in cross-check)
};

std::string formatMMSS(double seconds) {
    int s = static_cast<int>(seconds + 0.5);
    int m = s / 60;
    s = s % 60;
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << m << ":"
        << std::setfill('0') << std::setw(2) << s;
    return oss.str();
}

// Serialize one CrossCheckRow to the open CSV stream and flush so the line
// survives process kill. Keeps the same column layout as writeCrossCheckCSV
// for the OK/TIMEOUT/ERROR cases, minus DiffFromBaseline (computed offline).
void writeRowLine(std::ofstream& file, const CrossCheckRow& r) {
    std::string rel_str;
    if (r.status == CrossCheckStatus::OK) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(15) << r.reliability;
        rel_str = oss.str();
    }
    std::string time_str;
    if (r.status != CrossCheckStatus::TIMEOUT) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << r.time_seconds;
        time_str = oss.str();
    }
    std::string err = r.error_message;
    std::replace(err.begin(), err.end(), ',', ' ');
    std::replace(err.begin(), err.end(), '"', ' ');
    auto statusStr = [](CrossCheckStatus s) {
        switch (s) {
            case CrossCheckStatus::OK:      return "OK";
            case CrossCheckStatus::TIMEOUT: return "TIMEOUT";
            case CrossCheckStatus::ERROR:   return "ERROR";
        }
        return "UNKNOWN";
    };
    file << r.graph << ","
         << r.s << ","
         << r.t << ","
         << r.diameter << ","
         << r.method_id << ","
         << r.method_name << ","
         << statusStr(r.status) << ","
         << rel_str << ","
         << time_str << ","
         << r.recursions << ","
         << err << "\n";
    file.flush();
}

} // namespace

bool TestSuite::runCrossCheck(int timeout_sec,
                              const std::string& output_filename,
                              bool include_real_networks,
                              double tolerance,
                              int d_count,
                              int d_step,
                              const std::vector<int>& active_methods) {
    // timeout_sec is kept as a knob (CLI --timeout N) acting as an
    // override-all: if N > 0 and N != 30 (the default in main.cpp), the
    // caller wants uniform timeouts; otherwise we use kMethodTimeoutsSec.
    // We interpret "default 30" as "use tiered". Callers that want a
    // uniform override can pass anything != 30.
    const bool use_tiered = (timeout_sec == 30);
    LOG_INFO("Starting cross-check: timeout={}s (tiered={}), output={}, real_nets={}, tol={}, d_count={}, d_step={}",
             timeout_sec, use_tiered, output_filename, include_real_networks, tolerance,
             d_count, d_step);

    std::vector<std::string> method_names = getAvailableMethods();

    // Build the set of (graph, s, t, diameter) cells and load each graph
    // exactly once into an in-memory template. Algorithms mutate the graph
    // in place, so we copy from the template before each run instead of
    // re-parsing from disk.
    struct Cell {
        std::string graph;
        int s, t, d;
    };
    std::vector<Cell> cells;
    std::map<std::string, ReliabilityGraph> graph_templates;

    auto ensure_template = [&](const std::string& filename) -> ReliabilityGraph* {
        auto it = graph_templates.find(filename);
        if (it != graph_templates.end()) return &it->second;
        auto loaded = data_importer_.loadKAOGraph(filename);
        if (!loaded) return nullptr;
        auto [ins, _] = graph_templates.emplace(filename, std::move(*loaded));
        return &ins->second;
    };

    // Auto-ranging: for each graph we measure dist(s,t) and test
    // { dist, dist+step, ..., dist+(count-1)*step } clipped to |V|-1. This
    // avoids degenerate d < dist cells where every method trivially returns
    // R = 0 ("tautological agreement" instead of a real cross-check).
    auto add_auto_ranging = [&](const std::string& filename) {
        auto* tmpl = ensure_template(filename);
        if (!tmpl) {
            std::cerr << "  Skipping " << filename << ": could not load\n";
            return;
        }
        auto [s, t] = findSourceAndTargetVertices(*tmpl);
        if (s < 0 || t < 0) {
            std::cerr << "  Skipping " << filename << ": no s/t found\n";
            return;
        }
        auto ds = chooseDiameters(*tmpl, s, t, d_count, d_step);
        if (ds.empty()) {
            std::cerr << "  Skipping " << filename << ": s=" << s << ", t=" << t
                      << " unreachable within |V|-1 hops\n";
            return;
        }
        for (int d : ds) cells.push_back({filename, s, t, d});
    };

    for (const auto& [filename, diameters] : test_configurations_) {
        (void)diameters;  // diameters list ignored; d is auto-ranged from dist(s,t)
        add_auto_ranging(filename);
    }

    if (include_real_networks) {
        add_auto_ranging("Geant2004_kao.txt");
        add_auto_ranging("IEEE-118-node_kao.txt");
    }

    // Print the cell plan so it is immediately visible what diameters are
    // tested and on which (s,t) pairs. Any unreachable / skipped graphs are
    // already logged above.
    std::cout << "Cell plan (d_count=" << d_count << ", d_step=" << d_step << "):\n";
    {
        std::string prev_graph;
        int prev_s = -1, prev_t = -1;
        std::vector<int> ds_for_group;
        auto flush_group = [&]() {
            if (prev_graph.empty() || ds_for_group.empty()) return;
            std::cout << "  " << prev_graph
                      << " s=" << prev_s << " t=" << prev_t
                      << " dist=" << ds_for_group.front()
                      << " → d ∈ {";
            for (size_t i = 0; i < ds_for_group.size(); ++i) {
                if (i) std::cout << ", ";
                std::cout << ds_for_group[i];
            }
            std::cout << "}\n";
            ds_for_group.clear();
        };
        for (const auto& c : cells) {
            if (c.graph != prev_graph || c.s != prev_s || c.t != prev_t) {
                flush_group();
                prev_graph = c.graph;
                prev_s = c.s;
                prev_t = c.t;
            }
            ds_for_group.push_back(c.d);
        }
        flush_group();
    }
    std::cout << "Total cells: " << cells.size() << "\n";

    // Group cells by (graph, s, t) so methods with a multi-diameter (CDF) API
    // can amortize one factorization over every diameter in the group instead
    // of re-factoring per d. m3 / m4 / m5 internally already build the full
    // R(d) vector and used to discard everything except results[d]; the
    // grouped path now keeps the whole vector and writes one row per d from
    // a single call.
    struct Group {
        std::string graph;
        int s, t;
        std::vector<int> d_values;  // sorted ascending
    };
    std::vector<Group> groups;
    {
        std::map<std::tuple<std::string,int,int>, size_t> group_idx;
        for (const auto& c : cells) {
            auto key = std::make_tuple(c.graph, c.s, c.t);
            auto it = group_idx.find(key);
            if (it == group_idx.end()) {
                group_idx[key] = groups.size();
                groups.push_back({c.graph, c.s, c.t, {c.d}});
            } else {
                groups[it->second].d_values.push_back(c.d);
            }
        }
        for (auto& g : groups)
            std::sort(g.d_values.begin(), g.d_values.end());
    }

    auto methods_str = [&]() {
        std::ostringstream oss;
        for (size_t i = 0; i < active_methods.size(); ++i) {
            if (i) oss << ",";
            oss << "m" << active_methods[i];
        }
        return oss.str();
    }();
    if (use_tiered) {
        std::cout << "Cross-check: " << groups.size() << " (graph,s,t) groups × "
                  << active_methods.size() << " methods (" << methods_str << ")"
                  << " (tiered timeouts m0=" << kMethodTimeoutsSec[0]
                  << " m1=" << kMethodTimeoutsSec[1]
                  << " m2=" << kMethodTimeoutsSec[2]
                  << " m3=" << kMethodTimeoutsSec[3]
                  << " m4=" << kMethodTimeoutsSec[4]
                  << " m5=" << kMethodTimeoutsSec[5] << "s)\n";
    } else {
        std::cout << "Cross-check: " << groups.size() << " (graph,s,t) groups × "
                  << active_methods.size() << " methods (" << methods_str << ")"
                  << " (uniform timeout " << timeout_sec << "s)\n";
    }
    std::cout << "Methods 3/4/5 use one factorization per group; "
                 "0/1/2 still call once per d.\n";

    // Incremental CSV: open once, write header, flush after each row. Any
    // progress survives a kill.
    std::ofstream csv(output_filename);
    if (!csv) {
        std::cerr << "Could not open output file: " << output_filename << "\n";
        return false;
    }
    csv << "Graph,S,T,D,MethodId,Method,Status,Reliability,TimeSec,Recursions,Error\n";
    csv.flush();

    std::vector<CrossCheckRow> rows;
    rows.reserve(cells.size() * active_methods.size());

    auto run_start = std::chrono::steady_clock::now();

    auto build_row = [&](const std::string& graph, int s, int t, int d,
                         int method_id) {
        CrossCheckRow row;
        row.graph       = graph;
        row.s           = s;
        row.t           = t;
        row.diameter    = d;
        row.method_id   = method_id;
        row.method_name = method_names[method_id];
        return row;
    };

    int group_idx = 0;
    for (const auto& g : groups) {
        ++group_idx;
        std::cout << "[" << group_idx << "/" << groups.size() << "] "
                  << g.graph << " s=" << g.s << " t=" << g.t << " d ∈ {";
        for (size_t i = 0; i < g.d_values.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << g.d_values[i];
        }
        std::cout << "}\n";

        const int d_max = g.d_values.back();

        for (int method_id : active_methods) {
            const int cell_timeout = use_tiered
                ? kMethodTimeoutsSec[method_id]
                : timeout_sec;

            const bool use_cdf = (method_id == 3 || method_id == 4 || method_id == 5)
                                 && g.d_values.size() >= 2;

            ReliabilityGraph* tmpl = &graph_templates.at(g.graph);

            if (use_cdf) {
                auto graph = std::make_unique<ReliabilityGraph>(*tmpl);
                auto* gp = graph.get();
                int s = g.s, t = g.t;
                std::function<ReliabilityCdfResult()> calc;
                switch (method_id) {
                    case 3: calc = [gp, s, t, d_max]() { return gp->calculateReliabilityCdfMDecomposition(s, t, d_max); }; break;
                    case 4: calc = [gp, s, t, d_max]() { return gp->calculateReliabilityCdfCancelaPetingi(s, t, d_max); }; break;
                    case 5: calc = [gp, s, t, d_max]() { return gp->calculateReliabilityCdfMDecompositionCPFM(s, t, d_max); }; break;
                }
                auto graph_keeper = std::make_shared<std::unique_ptr<ReliabilityGraph>>(std::move(graph));
                std::function<ReliabilityCdfResult()> calc_wrapped =
                    [calc, graph_keeper]() { return calc(); };

                std::string err;
                auto t0 = std::chrono::steady_clock::now();
                auto opt = runWithTimeoutCdf(calc_wrapped, cell_timeout, &err);
                auto t1 = std::chrono::steady_clock::now();

                if (opt.has_value()) {
                    double total_time = opt->execution_time_sec;
                    long long total_recs = opt->recursions;
                    std::cout << "    m" << method_id << " [multi-d, "
                              << std::fixed << std::setprecision(3)
                              << total_time << "s]: ";
                    for (size_t i = 0; i < g.d_values.size(); ++i) {
                        int d = g.d_values[i];
                        CrossCheckRow row = build_row(g.graph, g.s, g.t, d, method_id);
                        row.status       = CrossCheckStatus::OK;
                        row.reliability  = (d < (int)opt->cdf.size()) ? opt->cdf[d] : 0.0;
                        row.time_seconds = total_time;
                        row.recursions   = (i == 0) ? total_recs : 0;
                        if (i) std::cout << ", ";
                        std::cout << "d=" << d << " R=" << std::setprecision(12)
                                  << row.reliability;
                        writeRowLine(csv, row);
                        rows.push_back(row);
                    }
                    std::cout << "\n";
                } else {
                    bool is_error = !err.empty();
                    double dt = std::chrono::duration<double>(t1 - t0).count();
                    std::cout << "    m" << method_id << " [multi-d]: "
                              << (is_error ? "ERROR (" : "TIMEOUT (>")
                              << (is_error ? err : (std::to_string(cell_timeout) + "s"))
                              << ")\n";
                    for (int d : g.d_values) {
                        CrossCheckRow row = build_row(g.graph, g.s, g.t, d, method_id);
                        row.status = is_error ? CrossCheckStatus::ERROR
                                              : CrossCheckStatus::TIMEOUT;
                        if (is_error) {
                            row.error_message = err;
                            row.time_seconds = dt;
                        }
                        writeRowLine(csv, row);
                        rows.push_back(row);
                    }
                }
            } else {
                // Per-d single calls (m0/m1/m2 always; m3/m4/m5 with single d).
                for (int d : g.d_values) {
                    auto graph = std::make_unique<ReliabilityGraph>(*tmpl);
                    auto* gp = graph.get();
                    int s = g.s, t = g.t;
                    std::function<ReliabilityResult()> calc;
                    switch (method_id) {
                        case 0: calc = [gp, s, t, d]() { return gp->calculateReliabilityBetweenVertices(s, t, d); }; break;
                        case 1: calc = [gp, s, t, d]() { return gp->calculateReliabilityWithRecursiveDecomposition(s, t, d); }; break;
                        case 2: calc = [gp, s, t, d]() { return gp->calculateReliabilityWithDecomposition(s, t, d); }; break;
                        case 3: calc = [gp, s, t, d]() { return gp->calculateReliabilityWithMDecomposition(s, t, d); }; break;
                        case 4: calc = [gp, s, t, d]() { return gp->calculateReliabilityCancelaPetingi(s, t, d); }; break;
                        case 5: calc = [gp, s, t, d]() { return gp->calculateReliabilityWithMDecompositionCPFM(s, t, d); }; break;
                    }
                    auto graph_keeper = std::make_shared<std::unique_ptr<ReliabilityGraph>>(std::move(graph));
                    std::function<ReliabilityResult()> calc_wrapped =
                        [calc, graph_keeper]() { return calc(); };

                    std::string err;
                    auto t0 = std::chrono::steady_clock::now();
                    auto opt = runWithTimeout(calc_wrapped, cell_timeout, &err);
                    auto t1 = std::chrono::steady_clock::now();

                    CrossCheckRow row = build_row(g.graph, g.s, g.t, d, method_id);
                    if (opt.has_value()) {
                        row.status       = CrossCheckStatus::OK;
                        row.reliability  = opt->reliability;
                        row.time_seconds = opt->execution_time_sec;
                        row.recursions   = opt->recursions;
                        std::cout << "    m" << method_id << " d=" << d
                                  << ": R=" << std::fixed << std::setprecision(12)
                                  << row.reliability
                                  << " (" << std::setprecision(3)
                                  << row.time_seconds << "s)\n";
                    } else if (!err.empty()) {
                        row.status = CrossCheckStatus::ERROR;
                        row.error_message = err;
                        row.time_seconds = std::chrono::duration<double>(t1 - t0).count();
                        std::cout << "    m" << method_id << " d=" << d
                                  << ": ERROR (" << err << ")\n";
                    } else {
                        row.status = CrossCheckStatus::TIMEOUT;
                        std::cout << "    m" << method_id << " d=" << d
                                  << ": TIMEOUT (>" << cell_timeout << "s)\n";
                    }
                    writeRowLine(csv, row);
                    rows.push_back(row);
                }
            }
        }

        // Progress / ETA after finishing all methods on this group.
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - run_start).count();
        double avg = elapsed / group_idx;
        double remaining = avg * (groups.size() - group_idx);
        std::cout << "    [" << group_idx << "/" << groups.size()
                  << " groups, elapsed=" << formatMMSS(elapsed)
                  << ", avg=" << std::fixed << std::setprecision(1) << avg
                  << "s/group, ETA=" << formatMMSS(remaining) << "]\n";
    }

    // Summary.
    int ok = 0, to = 0, er = 0;
    for (const auto& r : rows) {
        if (r.status == CrossCheckStatus::OK) ++ok;
        else if (r.status == CrossCheckStatus::TIMEOUT) ++to;
        else ++er;
    }
    std::cout << "\nSummary: total=" << rows.size()
              << " OK=" << ok << " TIMEOUT=" << to << " ERROR=" << er << "\n";

    // Find discrepancies (online in-memory baseline; the CSV already has all
    // rows without DiffFromBaseline — that column is computed in postprocess).
    auto key = [](const CrossCheckRow& r) {
        return r.graph + "|" + std::to_string(r.diameter);
    };
    std::map<std::string, double> baseline;
    for (int preferred : {0, 3, 5}) {
        for (const auto& r : rows) {
            if (r.method_id != preferred) continue;
            if (r.status != CrossCheckStatus::OK) continue;
            baseline.emplace(key(r), r.reliability);
        }
    }
    for (const auto& r : rows) {
        if (r.status == CrossCheckStatus::OK)
            baseline.emplace(key(r), r.reliability);
    }

    int discrepancies = 0;
    for (const auto& r : rows) {
        if (r.status != CrossCheckStatus::OK) continue;
        auto it = baseline.find(key(r));
        if (it == baseline.end()) continue;
        double diff = std::abs(r.reliability - it->second);
        if (diff > tolerance) {
            ++discrepancies;
            std::cout << "  MISMATCH " << r.graph << " d=" << r.diameter
                      << " m" << r.method_id << ": R=" << std::fixed
                      << std::setprecision(15) << r.reliability
                      << " baseline=" << it->second
                      << " diff=" << std::scientific << std::setprecision(3)
                      << diff << "\n";
        }
    }

    if (discrepancies == 0)
        std::cout << "No discrepancies > " << tolerance << ".\n";
    else
        std::cout << discrepancies << " discrepancies exceed tolerance "
                  << tolerance << ".\n";

    csv.close();
    std::cout << "Cross-check results written to " << output_filename << "\n";
    return discrepancies == 0;
}

} // namespace graph_reliability

