#include "DataImporter.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cctype>

DataImporter::DataImporter(const std::string& baseDataPath) 
    : baseDataPath_(normalizePath(baseDataPath)) {
    // Убеждаемся, что путь заканчивается разделителем
    if (!baseDataPath_.empty() && baseDataPath_.back() != '/' && baseDataPath_.back() != '\\') {
        baseDataPath_ += '/';
    }
}

std::unique_ptr<kGraph> DataImporter::loadKAOGraph(const std::string& filename) {
    std::string fullPath = getFullPath(filename);
    
    if (!fileExists(filename)) {
        throw std::runtime_error("File not found: " + fullPath);
    }
    
    std::ifstream fin(fullPath);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + fullPath);
    }
    
    try {
        auto graph = loadKAOFromStream(fin);
        validateGraph(graph.get(), filename);
        return graph;
    } catch (const std::exception& e) {
        throw std::runtime_error("Error loading KAO graph from " + filename + ": " + e.what());
    }
}

std::unique_ptr<kGraph> DataImporter::loadEdgeListGraph(const std::string& filename, double reliability) {
    std::string fullPath = getFullPath(filename);
    
    if (!fileExists(filename)) {
        throw std::runtime_error("File not found: " + fullPath);
    }
    
    std::ifstream fin(fullPath);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + fullPath);
    }
    
    try {
        // Читаем рёбра
        std::vector<std::pair<int, int>> edges;
        int u, v;
        char dash1, dash2;
        
        while (fin >> u >> dash1 >> dash2 >> v) {
            if (dash1 != '-' || dash2 != '-') {
                throw std::runtime_error("Invalid edge format in file: " + filename);
            }
            edges.emplace_back(u, v);
        }
        
        if (edges.empty()) {
            throw std::runtime_error("No edges found in file: " + filename);
        }
        
        // Определяем количество вершин
        int maxVertex = 0;
        for (const auto& edge : edges) {
            maxVertex = std::max(maxVertex, std::max(edge.first, edge.second));
        }
        int n = maxVertex;
        
        // Строим список смежности
        std::vector<std::vector<int>> adj(n + 1);
        for (const auto& edge : edges) {
            adj[edge.first].push_back(edge.second);
            adj[edge.second].push_back(edge.first);
        }
        
        // Сортируем списки смежности
        for (int i = 1; i <= n; ++i) {
            std::sort(adj[i].begin(), adj[i].end());
        }
        
        // Строим CSR формат (KAO, FO)
        std::vector<int> KAO(n + 2, 0);
        std::vector<int> FO;
        
        for (int i = 1; i <= n; ++i) {
            KAO[i + 1] = KAO[i] + static_cast<int>(adj[i].size());
        }
        
        FO.reserve(KAO[n + 1]);
        for (int i = 1; i <= n; ++i) {
            for (int neighbor : adj[i]) {
                FO.push_back(neighbor);
            }
        }
        
        // Targets (все 0, кроме первой и последней вершины)
        std::vector<int> Targets(n + 1, 0);
        if (n >= 1) {
            Targets[1] = 1;  // Первая вершина - источник
            Targets[n] = 1;  // Последняя вершина - приёмник
        }
        
        // PArray (все рёбра имеют одинаковую надёжность)
        std::vector<double> PArray(FO.size(), reliability);
        
        auto graph = std::make_unique<kGraph>(KAO, FO, PArray, Targets);
        validateGraph(graph.get(), filename);
        
        return graph;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error loading Edge List graph from " + filename + ": " + e.what());
    }
}

std::vector<std::unique_ptr<kGraph>> DataImporter::loadGraphsToMerge() {
    std::string filename = "GraphsToMerge.txt";
    std::string fullPath = getFullPath(filename);
    
    if (!fileExists(filename)) {
        throw std::runtime_error("File not found: " + fullPath);
    }
    
    std::ifstream fin(fullPath);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + fullPath);
    }
    
    std::vector<std::unique_ptr<kGraph>> graphs;
    
    try {
        // Читаем графы до конца файла
        while (fin.good() && !fin.eof()) {
            // Пропускаем пустые строки и комментарии
            std::string line;
            std::getline(fin, line);
            
            // Убираем пробелы в начале и конце
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            // Пропускаем пустые строки и комментарии
            if (line.empty() || line[0] == '/' || line[0] == '#') {
                continue;
            }
            
            // Возвращаемся к началу строки для чтения графа
            fin.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            
            auto graph = loadKAOFromStream(fin);
            validateGraph(graph.get(), filename);
            graphs.push_back(std::move(graph));
        }
        
        if (graphs.empty()) {
            throw std::runtime_error("No valid graphs found in " + filename);
        }
        
        return graphs;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error loading graphs from " + filename + ": " + e.what());
    }
}

std::string DataImporter::getFullPath(const std::string& filename) const {
    return baseDataPath_ + filename;
}

bool DataImporter::fileExists(const std::string& filename) const {
    std::string fullPath = getFullPath(filename);
    std::ifstream file(fullPath);
    return file.good();
}

std::vector<std::string> DataImporter::getAvailableFiles() const {
    std::vector<std::string> files;
    
    try {
        // Используем std::filesystem если доступен (C++17)
        #if __cplusplus >= 201703L
        for (const auto& entry : std::filesystem::directory_iterator(baseDataPath_)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                // Фильтруем только текстовые файлы
                if (filename.length() >= 4 && 
                    (filename.substr(filename.length() - 4) == ".txt" ||
                     filename.substr(filename.length() - 4) == ".kao")) {
                    files.push_back(filename);
                }
            }
        }
        #else
        // Fallback для старых версий C++
        // Здесь можно добавить альтернативную реализацию
        #endif
    } catch (const std::exception& e) {
        std::cerr << "Warning: Cannot list directory " << baseDataPath_ 
                  << ": " << e.what() << std::endl;
    }
    
    std::sort(files.begin(), files.end());
    return files;
}

void DataImporter::setBasePath(const std::string& newBasePath) {
    baseDataPath_ = normalizePath(newBasePath);
    if (!baseDataPath_.empty() && baseDataPath_.back() != '/' && baseDataPath_.back() != '\\') {
        baseDataPath_ += '/';
    }
}

void DataImporter::convertEdgeListToKAO(const std::string& inputFilename, 
                                        const std::string& outputFilename, 
                                        double reliability) {
    try {
        auto graph = loadEdgeListGraph(inputFilename, reliability);
        
        std::string outputPath = getFullPath(outputFilename);
        std::ofstream fout(outputPath);
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot create output file: " + outputPath);
        }
        
        // Записываем в KAO формате
        graph->kGraphFileOutput(fout);
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error converting Edge List to KAO: " + std::string(e.what()));
    }
}

void DataImporter::convertKAOToEdgeList(const std::string& inputFilename, 
                                        const std::string& outputFilename) {
    try {
        auto graph = loadKAOGraph(inputFilename);
        
        std::string outputPath = getFullPath(outputFilename);
        std::ofstream fout(outputPath);
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot create output file: " + outputPath);
        }
        
        // Конвертируем в Edge List формат
        graph->convertKAOFOToEdgeList(inputFilename, outputPath);
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error converting KAO to Edge List: " + std::string(e.what()));
    }
}

std::unique_ptr<kGraph> DataImporter::loadKAOFromStream(std::ifstream& fin) {
    if (!fin.good()) {
        throw std::runtime_error("Input stream is not in good state");
    }
    
    std::string line;
    
    // Читаем KAO
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Cannot read KAO line");
    }
    
    std::istringstream skao(line);
    std::string temp;
    std::vector<int> KAO;
    
    while (std::getline(skao, temp, ',')) {
        try {
            KAO.push_back(std::stoi(temp));
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid KAO value: " + temp);
        }
    }
    
    if (KAO.empty()) {
        throw std::runtime_error("KAO array is empty");
    }
    
    // Читаем FO
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Cannot read FO line");
    }
    
    std::istringstream sfo(line);
    std::vector<int> FO;
    
    while (std::getline(sfo, temp, ',')) {
        try {
            FO.push_back(std::stoi(temp));
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid FO value: " + temp);
        }
    }
    
    // Читаем Targets
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Cannot read Targets line");
    }
    
    std::istringstream stargets(line);
    std::vector<int> Targets;
    
    while (std::getline(stargets, temp, ',')) {
        try {
            Targets.push_back(std::stoi(temp));
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid Target value: " + temp);
        }
    }
    
    // Читаем PArray (надежность)
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Cannot read PArray line");
    }
    
    // Убираем комментарии
    size_t commentPos = line.find("//");
    if (commentPos != std::string::npos) {
        line = line.substr(0, commentPos);
    }
    
    // Убираем пробелы
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    
    if (line.empty()) {
        throw std::runtime_error("PArray line is empty");
    }
    
    // Заменяем запятую на точку для корректного парсинга
    size_t found = line.find(',');
    if (found != std::string::npos) {
        line[found] = '.';
    }
    
    double p;
    try {
        p = std::stod(line);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid PArray value: " + line);
    }
    
    if (p < 0.0 || p > 1.0) {
        throw std::runtime_error("PArray value out of range [0,1]: " + std::to_string(p));
    }
    
    std::vector<double> PArray(FO.size(), p);
    
    return std::make_unique<kGraph>(KAO, FO, PArray, Targets);
}

void DataImporter::validateGraph(const kGraph* graph, const std::string& filename) const {
    if (!graph) {
        throw std::runtime_error("Graph is null");
    }
    
    // Проверяем размеры массивов
    if (graph->KAO.empty()) {
        throw std::runtime_error("KAO array is empty in " + filename);
    }
    
    if (graph->FO.empty()) {
        throw std::runtime_error("FO array is empty in " + filename);
    }
    
    if (graph->PArray.empty()) {
        throw std::runtime_error("PArray is empty in " + filename);
    }
    
    if (graph->Targets.empty()) {
        throw std::runtime_error("Targets array is empty in " + filename);
    }
    
    // Проверяем соответствие размеров
    if (graph->FO.size() != graph->PArray.size()) {
        throw std::runtime_error("FO and PArray sizes don't match in " + filename);
    }
    
    if (graph->KAO.size() != graph->Targets.size()) {
        throw std::runtime_error("KAO and Targets sizes don't match in " + filename);
    }
    
    // Проверяем корректность KAO
    if (graph->KAO[0] != 0) {
        throw std::runtime_error("KAO[0] must be 0 in " + filename);
    }
    
    if (graph->KAO.back() != static_cast<int>(graph->FO.size())) {
        throw std::runtime_error("KAO last element doesn't match FO size in " + filename);
    }
    
    // Проверяем монотонность KAO
    for (size_t i = 1; i < graph->KAO.size(); ++i) {
        if (graph->KAO[i] < graph->KAO[i-1]) {
            throw std::runtime_error("KAO array is not monotonic in " + filename);
        }
    }
    
    // Проверяем диапазон значений в FO
    int numVertices = static_cast<int>(graph->KAO.size()) - 1;
    for (int vertex : graph->FO) {
        if (vertex < 1 || vertex > numVertices) {
            throw std::runtime_error("Invalid vertex index in FO array in " + filename);
        }
    }
    
    // Проверяем диапазон значений в PArray
    for (double prob : graph->PArray) {
        if (prob < 0.0 || prob > 1.0) {
            throw std::runtime_error("Invalid probability value in PArray in " + filename);
        }
    }
}

std::string DataImporter::normalizePath(const std::string& path) const {
    std::string normalized = path;
    
    // Заменяем обратные слеши на прямые (для кроссплатформенности)
    std::replace(normalized.begin(), normalized.end(), '\\', '/');
    
    // Убираем лишние пробелы
    normalized.erase(0, normalized.find_first_not_of(" \t\r\n"));
    normalized.erase(normalized.find_last_not_of(" \t\r\n") + 1);
    
    return normalized;
}
