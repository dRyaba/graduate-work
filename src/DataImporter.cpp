/**
 * @file DataImporter.cpp
 * @brief Implementation of the DataImporter class
 * @author Graduate Work Project
 * @date 2026
 */

#include "graph_reliability/DataImporter.h"
#include "graph_reliability/Logger.h"
#include "graph_reliability/Exceptions.h"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace graph_reliability {

DataImporter::DataImporter(const std::string& base_data_path) 
    : base_data_path_(base_data_path) {
    // Ensure base path ends with separator
    if (!base_data_path_.empty() && base_data_path_.back() != '/' && base_data_path_.back() != '\\') {
        base_data_path_ += '/';
    }
    LOG_DEBUG("DataImporter initialized with base path: {}", base_data_path_);
}

std::string DataImporter::getFullPath(const std::string& filename) const {
    // If filename is already absolute, return it
    fs::path p(filename);
    if (p.is_absolute()) {
        return filename;
    }
    return base_data_path_ + filename;
}

bool DataImporter::fileExists(const std::string& filename) const {
    return fs::exists(getFullPath(filename));
}

std::unique_ptr<ReliabilityGraph> DataImporter::loadKAOGraph(const std::string& filename) {
    std::string path = getFullPath(filename);
    LOG_DEBUG("Loading KAO graph from: {}", path);
    std::ifstream file(path);
    if (!file) {
        LOG_ERROR("Cannot open KAO graph file: {}", path);
        throw FileNotFoundException(path);
    }
    auto graph = loadKAOFromStream(file);
    LOG_DEBUG("Successfully loaded KAO graph: {} vertices, {} edges", 
              graph->numVertices(), graph->numEdges());
    return graph;
}

std::unique_ptr<ReliabilityGraph> DataImporter::loadKAOFromStream(std::ifstream& input_stream) {
    LOG_DEBUG("Parsing KAO format from stream");
    std::string line;
    
    // 1. Read KAO (Offsets)
    if (!std::getline(input_stream, line)) {
        LOG_ERROR("Empty file or missing KAO line");
        throw InvalidFormatException("Empty file or missing KAO line");
    }
    std::vector<Graph::VertexId> kao;
    std::stringstream ss_kao(line);
    std::string token;
    while (std::getline(ss_kao, token, ',')) {
        kao.push_back(std::stoi(token));
    }

    // 2. Read FO (Adjacency)
    if (!std::getline(input_stream, line)) {
        LOG_ERROR("Missing FO line");
        throw InvalidFormatException("Missing FO line");
    }
    std::vector<Graph::VertexId> fo;
    std::stringstream ss_fo(line);
    while (std::getline(ss_fo, token, ',')) {
        fo.push_back(std::stoi(token));
    }

    // 3. Read Targets
    if (!std::getline(input_stream, line)) {
        LOG_ERROR("Missing Targets line");
        throw InvalidFormatException("Missing Targets line");
    }
    std::vector<int> targets;
    std::stringstream ss_targets(line);
    while (std::getline(ss_targets, token, ',')) {
        targets.push_back(std::stoi(token));
    }

    // 4. Read Probability
    // Legacy format: "0.9" (single value) or list?
    // Legacy kGraphFileInput reads one line. 
    // "size_t found = line.find(','); if (found != npos) ... double p = stod(line);"
    // It creates Parray of size FO.size() filled with p.
    // OR it might support per-edge P?
    // Legacy code: `std::vector<double> Parray(FO.size(), p);`
    // So it expects a single value for all edges usually.
    
    if (!std::getline(input_stream, line)) {
        LOG_ERROR("Missing Probability line");
        throw InvalidFormatException("Missing Probability line");
    }
    
    // Handle comma vs dot? Legacy does `line[found] = '.'`.
    size_t found = line.find(',');
    if (found != std::string::npos) {
        line[found] = '.';
    }
    
    double reliability = std::stod(line);
    LOG_DEBUG("Loaded reliability value: {}", reliability);
    std::vector<Graph::Probability> p_array(fo.size(), reliability);

    // 5. Create ReliabilityGraph
    // Note: Legacy data is 1-based? 
    // If my ReliabilityGraph implementation assumes 0-based, I should convert here?
    // Yes, for Professional API, we should normalize data at the boundary (Importer).
    
    // Convert KAO: KAO in file starts with 0 (if v1 starts at 0).
    // If file KAO is 0, 2, ... 
    // My 0-based Graph expects KAO[0]=0.
    // So KAO structure is compatible if I map v -> v-1.
    
    // Convert FO: File has 1-based indices (e.g., 2).
    // I need to convert to 0-based (e.g., 1).
    for (auto& v : fo) {
        if (v > 0) v--; // Decrement to 0-based
    }
    
    // Convert Targets: File has vector of 0/1. Size N+1 or N?
    // Legacy Targets is size N+1 (index 0 unused).
    // My 0-based Graph uses 0..N-1.
    // So I should remove first element?
    if (!targets.empty()) {
        targets.erase(targets.begin());
    }

    return std::make_unique<ReliabilityGraph>(std::move(kao), std::move(fo), std::move(p_array), std::move(targets));
}

std::unique_ptr<ReliabilityGraph> DataImporter::loadEdgeListGraph(const std::string& filename, double reliability) {
    // Convert EdgeList to temporary KAO logic then load?
    // Or just reuse Graph::convert... but that writes to file.
    // I should implement direct loading.
    // For now, let's use the file conversion utility as a bridge if needed, 
    // or just parse edge list directly.
    
    std::string path = getFullPath(filename);
    LOG_DEBUG("Loading EdgeList graph from: {} with reliability={}", path, reliability);
    std::ifstream file(path);
    if (!file) {
        LOG_ERROR("Cannot open EdgeList file: {}", path);
        throw FileNotFoundException(path);
    }
    
    // Basic parsing logic similar to Graph::convertEdgeListToKAOFO
    std::vector<std::pair<int, int>> edges;
    std::string line;
    int u, v;
    int max_v = 0;
    
    while (std::getline(file, line)) {
        size_t comment = line.find("//");
        if (comment != std::string::npos) line = line.substr(0, comment);
        if (line.empty()) continue;
        
        for(char& c : line) if(!isdigit(c)) c = ' ';
        std::stringstream ss(line);
        if (ss >> u >> v) {
            if (u > 0 && v > 0) { // Assume 1-based input
                edges.emplace_back(u - 1, v - 1); // Convert to 0-based
                max_v = std::max(max_v, std::max(u-1, v-1));
            }
        }
    }
    
    int n = max_v + 1;
    std::vector<std::vector<int>> adj(n);
    for (const auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    
    std::vector<Graph::VertexId> kao(n + 1);
    kao[0] = 0;
    for (int i = 0; i < n; ++i) {
        std::sort(adj[i].begin(), adj[i].end());
        kao[i + 1] = kao[i] + adj[i].size();
    }
    
    std::vector<Graph::VertexId> fo;
    fo.reserve(kao[n]);
    for (const auto& neighbors : adj) {
        for (int neighbor : neighbors) {
            fo.push_back(neighbor);
        }
    }
    
    std::vector<Graph::Probability> p_array(fo.size(), reliability);
    std::vector<int> targets(n, 0); // Default targets? Or none?
    
    LOG_DEBUG("Successfully loaded EdgeList graph: {} vertices, {} edges", n, edges.size() * 2);
    return std::make_unique<ReliabilityGraph>(std::move(kao), std::move(fo), std::move(p_array), std::move(targets));
}

std::vector<std::unique_ptr<ReliabilityGraph>> DataImporter::loadGraphsToMerge() {
    std::string list_file = getFullPath("GraphsToMerge.txt");
    LOG_INFO("Loading graphs from merge list: {}", list_file);
    std::ifstream file(list_file);
    if (!file) {
        LOG_ERROR("Cannot open GraphsToMerge.txt: {}", list_file);
        throw FileNotFoundException(list_file);
    }
    
    std::vector<std::unique_ptr<ReliabilityGraph>> graphs;
    std::string filename;
    while (std::getline(file, filename)) {
        if (!filename.empty()) {
             // Clean filename (remove \r etc)
             filename.erase((std::remove)(filename.begin(), filename.end(), '\r'), filename.end());
             LOG_DEBUG("Loading graph for merge: {}", filename);
             graphs.push_back(loadKAOGraph(filename));
        }
    }
    LOG_INFO("Loaded {} graphs for merging", graphs.size());
    return graphs;
}

std::vector<std::string> DataImporter::getAvailableFiles() const {
    std::vector<std::string> files;
    if (fs::exists(base_data_path_)) {
        for (const auto& entry : fs::directory_iterator(base_data_path_)) {
            if (entry.path().extension() == ".txt" || entry.path().extension() == ".kao") {
                files.push_back(entry.path().filename().string());
            }
        }
    }
    return files;
}

void DataImporter::setBasePath(const std::string& new_base_path) {
    base_data_path_ = new_base_path;
    if (!base_data_path_.empty() && base_data_path_.back() != '/' && base_data_path_.back() != '\\') {
        base_data_path_ += '/';
    }
}

void DataImporter::convertEdgeListToKAO(const std::string& input_filename, 
                                       const std::string& output_filename, 
                                       double reliability) {
    LOG_INFO("Converting EdgeList to KAO: {} -> {} (reliability={})", 
             input_filename, output_filename, reliability);
    // Use temporary Graph to convert
    // Since convertEdgeListToKAOFO is member of Graph and we need to use it:
    Graph g; // Dummy
    g.convertEdgeListToKAOFO(getFullPath(input_filename), getFullPath(output_filename), reliability);
    LOG_INFO("Conversion completed successfully");
}

void DataImporter::convertKAOToEdgeList(const std::string& input_filename, 
                                       const std::string& output_filename) {
    LOG_INFO("Converting KAO to EdgeList: {} -> {}", input_filename, output_filename);
    Graph g;
    g.convertKAOFOToEdgeList(getFullPath(input_filename), getFullPath(output_filename));
    LOG_INFO("Conversion completed successfully");
}

void DataImporter::validateGraph(const ReliabilityGraph* graph, const std::string& filename) const {
    if (!graph) {
        LOG_ERROR("Graph validation failed: graph is null for file {}", filename);
        throw InvalidGraphException("Graph is null: " + filename);
    }
    if (graph->numVertices() == 0 && graph->numEdges() > 0) {
        LOG_ERROR("Graph validation failed: edges without vertices for file {}", filename);
        throw InvalidGraphException("Edges without vertices: " + filename);
    }
    LOG_DEBUG("Graph validation passed for file {}", filename);
}

std::string DataImporter::normalizePath(const std::string& path) const {
    std::string res = path;
    std::replace(res.begin(), res.end(), '\\', '/');
    return res;
}

} // namespace graph_reliability

