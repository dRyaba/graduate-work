#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <chrono>
#include <numeric>
#include <algorithm>
#include "kGraphOperations.h"
#include "DataImporter.h"
#include <map>
#include <iomanip>
#include <string>


int NUMBER_OF_REPETITIONS = 1;

// —— структуры для тестов ——————————————————————————————

struct TestConfig {
    std::string graph_file_name;
    int s_node;               // 1‑based номер источника
    int t_node;               // 1‑based номер приёмника
    int upper_bound_diameter;
    int num_repetitions;
};

struct TestRunResult {
    double time_sec;
    double reliability;
    long long factoring_recursions;
};

// —— функции конвертации (используют DataImporter) ——————————————————————————————

// Из edge‑list (строка "u -- v") → KAO/FO+Targets+PArray
void convertEdgeListToKAOFO(const std::string& inPath,
    const std::string& outPath,
    double reliability)
{
    try {
        DataImporter importer;
        importer.convertEdgeListToKAO(inPath, outPath, reliability);
        std::cout << "EdgeList -> KAOFO written to " << outPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error converting EdgeList to KAO: " << e.what() << "\n";
    }
}

// Из KAO/FO(+…) → edge‑list ("u -- v")
void convertKAOFOToEdgeList(const std::string& inPath,
    const std::string& outPath)
{
    try {
        DataImporter importer;
        importer.convertKAOToEdgeList(inPath, outPath);
        std::cout << "KAOFO → EdgeList written to " << outPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error converting KAO to EdgeList: " << e.what() << "\n";
    }
}

// —— тестовая функция  ——————————————————————————

TestRunResult run_single_test_config(const TestConfig& cfg,
    int methodID /* 0=MDecomp,1=Rec,2=Simple */)
{
    DataImporter importer("graphs_data/");
    std::vector<double> times;
    std::vector<double> reliabs;
    std::vector<long long> recs;

    for (int i = 0; i < cfg.num_repetitions; ++i) {
        try {
            auto graphPtr = importer.loadKAOGraph(cfg.graph_file_name);
            kGraph G = *graphPtr; // Копируем граф

            //G.CutPointsSearch(52, -1); //в переменной cutpoints будет лежать список точек сочленения

            auto t0 = std::chrono::high_resolution_clock::now();
            GraphMethodResult mr;
            if (methodID == 0) {
                mr = G.ReliabilityDiamConstr2VertMDecompose(
                    cfg.s_node, cfg.t_node, cfg.upper_bound_diameter);
            }
            else if (methodID == 1) {
                mr = G.ReliabilityDiamConstr2VertRecursiveDecomposition(
                    cfg.s_node, cfg.t_node, cfg.upper_bound_diameter);
            }
            else if (methodID == 2) {
                mr = G.ReliabilityDiamConstr2VertDecomposeSimpleFacto(
                    cfg.s_node, cfg.t_node, cfg.upper_bound_diameter);
            }
            else if (methodID == 3) {
                mr = G.ReliabilityDiamConstr2Vert(
                    cfg.s_node, cfg.t_node, cfg.upper_bound_diameter);
            }
            auto t1 = std::chrono::high_resolution_clock::now();

            times.push_back(
                std::chrono::duration<double>(t1 - t0).count()
            );
            reliabs.push_back(mr.reliability);
            recs.push_back(mr.recursions);
        } catch (const std::exception& e) {
            std::cerr << "Error loading graph " << cfg.graph_file_name 
                      << ": " << e.what() << std::endl;
            // Возвращаем нулевой результат при ошибке
            times.push_back(0.0);
            reliabs.push_back(0.0);
            recs.push_back(0);
        }
    }

    TestRunResult out{};
    out.time_sec = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    out.reliability = std::accumulate(reliabs.begin(), reliabs.end(), 0.0) / reliabs.size();
    out.factoring_recursions = static_cast<long long>(
        std::accumulate(recs.begin(), recs.end(), 0LL) / recs.size()
        );
    return out;
}

// —— main с режимами —————————————————————————————

int main(int argc, char* argv[]) {
    // Для настройки параметров запуска необходимо: 
    // Находясь в файле graduate.cpp Проект -> свойства -> Свойства конфигурации -> Отладка -> 
    // -> Аргументы команды: --convert edge2kao 4_blocks_sausage_3x3_edgelist.txt 4_blocks_sausage_3x3_kao.txt 0.9

    if (argc >= 2 && std::string(argv[1]) == "--convert") {
        if (argc < 5) {
            std::cerr << "Usage:\n  --convert edge2kao <in> <out> <reliab>\n"
                "  --convert kao2edge <in> <out>\n";
            return 1;
        }
        std::string dir = argv[2];
        std::string in = argv[3];
        std::string out = argv[4];
        if (dir == "edge2kao") {
            if (argc < 6) {
                std::cerr << "Missing reliability value\n";
                return 1;
            }
            double r = std::stod(argv[5]);
            convertEdgeListToKAOFO(in, out, r);
        }
        else if (dir == "kao2edge") {
            convertKAOFOToEdgeList(in, out);
        }
        else {
            std::cerr << "Unknown convert mode: " << dir << "\n";
            return 1;
        }
        return 0;
    }

    // —— иначе — запускаем серию тестов ————————————————
    std::cout << "Running tests...\n";
    
    // Инициализируем DataImporter
    DataImporter importer("graphs_data/");
    
    std::map<std::string, std::vector<int>> fileDiameters = {
        //{ "K4_kao.txt", { 1, 2, 3 } },
        //{ "2_3x3_blocks_kao.txt", { 8, 9, 10 } },
        { "3_blocks_sausage_3x3_kao.txt", { 9, 10, 11, 12 } },
        { "4_blocks_sausage_3x3_kao.txt", { 13, 14, 15, 16 } },
        { "5_blocks_sausage_3x3_kao.txt", { 17, 18, 19, 20 } },
        { "6_blocks_sausage_3x3_kao.txt", { 21, 22, 23, 24 } },
        //{ "3_blocks_sausage_4x4_kao.txt", { 11, 13/*, 15, 17, 19*/ } },
        //{ "Geant2004_kao.txt", { 6, 7, 8, 9, 10, 11} },
        //{ "IEEE-118-node_kao.txt", { 8, 9/*, 10, 11 */} },
        //{ "UPS_of_Russia_composed_with_colored_cut_vertices_kao.txt", { 18, 19/*, 20, 21, 22, 23, 24, 25*/ } }
        //{ "Geant2009_kao.txt", { 8,  12, 16 } }
    };

    // для тестов выбираем метод по флагу:
    // 0 — MDecompose, 1 — RecursiveDecompSF(Migov's Formula), 2 — DecompSF(Ryabinin-Migov's Formula), 3 - SF
    int methodID = 3;  // или считывать из аргументов

    std::vector<TestConfig> suite;
    for (const auto& pair : fileDiameters) {
        const std::string& fname = pair.first;
        const std::vector<int>& diamVec = pair.second;

        int s = -1, t = -1;
        {
            try {
                auto graphPtr = importer.loadKAOGraph(fname);
                const kGraph& G = *graphPtr;
                
                // Ищем source и target вершины
                for (int i = 1; i < static_cast<int>(G.Targets.size()); ++i) {
                    if (G.Targets[i] == 1) {
                        if (s < 0) s = i;
                        else if (t < 0) t = i;
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Warning: cannot load " << fname << ": " << e.what() << "\n";
                continue;
            }
        }

        if (s < 1 || t < 1) {
            std::cerr << "Warning: no valid s/t found in " << fname << "\n";
            continue;
        }

        for (int D : diamVec) {
            suite.push_back({ fname, s, t, D, NUMBER_OF_REPETITIONS });
        }
    }

    std::ofstream sumF("test_summary_results.csv");
    sumF << "Graph,S,T,D,Reps,Method,AvgTime,AvgRel,AvgRecs\n";

    //bool useRec = true; // Можно добавить цикл для тестирования обоих методов
    for (const auto& cfg : suite) {
        std::cout << "Test " << cfg.graph_file_name
            << " s=" << cfg.s_node
            << " t=" << cfg.t_node
            << " D<=" << cfg.upper_bound_diameter;

        auto R = run_single_test_config(cfg, methodID);
        std::vector<std::string> methodNames = {
            "ReliabilityDiamConstr2VertMDecompose",
            "ReliabilityDiamConstr2VertRecursiveDecomposition",
            "ReliabilityDiamConstr2VertDecomposeSimpleFacto",
            "ReliabilityDiamConstr2VertSimpleFacto"
        };
        std::cout << " Method=" << methodNames[methodID] << ","
            << " Time=" << R.time_sec << ","
            << " Reliab=" << R.reliability << ","
            << " recursions=" << R.factoring_recursions << "\n";

        sumF << cfg.graph_file_name << ","
            << cfg.s_node << ","
            << cfg.t_node << ","
            << cfg.upper_bound_diameter << ","
            << cfg.num_repetitions << ","
            << methodNames[methodID] << ","
            << R.time_sec << ","
            << R.reliability << ","
            << R.factoring_recursions << "\n";
    }
    sumF.close();

    std::cout << "All done. Summary in test_summary_results.csv\n";
    return 0;
}