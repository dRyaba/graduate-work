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

// —— функции конвертации ——————————————————————————————

// Из edge‑list (строка "u -- v") → KAO/FO+Targets+PArray
void convertEdgeListToKAOFO(const std::string& inPath,
    const std::string& outPath,
    double reliability)
{
    std::ifstream fin(inPath);
    if (!fin) {
        std::cerr << "Error opening edge‑list file: " << inPath << "\n";
        return;
    }

    // 1) читаем рёбра
    std::vector<std::pair<int, int>> edges;
    int u, v; char dash;
    while (fin >> u >> dash >> dash >> v) {
        edges.emplace_back(u, v);
    }
    fin.close();

    // 2) определяем n = max{u,v}
    int maxV = 0;
    for (auto& e : edges) {
        maxV = std::max({ maxV, e.first, e.second });
    }
    int n = maxV;

    // 3) строим список смежности
    std::vector<std::vector<int>> adj(n + 1);
    for (auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    for (int i = 1; i <= n; ++i) {
        std::sort(adj[i].begin(), adj[i].end());
    }

    // 4) CSR: KAO и FO
    std::vector<int> KAO(n + 2, 0), FO;
    for (int i = 1; i <= n; ++i) {
        KAO[i + 1] = KAO[i] + static_cast<int>(adj[i].size());
    }
    FO.reserve(KAO[n + 1]);
    for (int i = 1; i <= n; ++i) {
        for (int w : adj[i]) {
            FO.push_back(w);
        }
    }

    // 5) Targets (все 0) и PArray (все = reliability)
    std::vector<int> Targets(n + 1, 0);
    double PArray = reliability;

    // 6) записываем: KAO, FO, Targets, PArray
    std::ofstream fout(outPath);
    if (!fout) {
        std::cerr << "Error opening output file: " << outPath << "\n";
        return;
    }
    // KAO
    for (int i = 1; i <= n + 1; ++i) {
        fout << KAO[i] << (i < n + 1 ? ',' : '\n');
    }
    // FO
    for (size_t i = 0; i < FO.size(); ++i) {
        fout << FO[i] << (i + 1 < FO.size() ? ',' : '\n');
    }
    // Targets
    for (int i = 0; i < n + 1; ++i) {
        fout << Targets[i] << (i + 1 < n + 1 ? ',' : '\n');
    }
    // PArray
    //for (size_t i = 0; i < PArray.size(); ++i) {
    //    fout << PArray[i] << (i + 1 < PArray.size() ? ',' : '\n');
    //}
    fout << PArray << '\n';

    fout.close();

    std::cout << "EdgeList -> KAOFO written to " << outPath << "\n";
}

// Из KAO/FO(+…) → edge‑list ("u -- v")
void convertKAOFOToEdgeList(const std::string& inPath,
    const std::string& outPath)
{
    std::ifstream fin(inPath);
    if (!fin) {
        std::cerr << "Error opening KAOFO file: " << inPath << "\n";
        return;
    }

    std::string line, token;
    std::vector<int> KAO, FO;

    // 1) читаем первую строку → KAO
    if (std::getline(fin, line)) {
        std::istringstream ss(line);
        while (std::getline(ss, token, ',')) {
            KAO.push_back(std::stoi(token));
        }
    }
    // 2) читаем вторую строку → FO
    if (std::getline(fin, line)) {
        std::istringstream ss(line);
        while (std::getline(ss, token, ',')) {
            FO.push_back(std::stoi(token));
        }
    }
    // пропускаем Targets, PArray
    std::getline(fin, line);
    std::getline(fin, line);
    fin.close();

    int n = static_cast<int>(KAO.size()) - 1;
    std::set<std::pair<int, int>> seen;
    std::ofstream fout(outPath);
    if (!fout) {
        std::cerr << "Error opening edge‑list output: " << outPath << "\n";
        return;
    }

    for (int u = 1; u <= n; ++u) {
        int start = KAO[u - 1], end = KAO[u];
        for (int idx = start; idx < end; ++idx) {
            int v = FO[idx];
            auto e = std::minmax(u, v);
            if (seen.insert(e).second) {
                fout << e.first << " -- " << e.second << "\n";
            }
        }
    }
    fout.close();
    std::cout << "KAOFO → EdgeList written to " << outPath << "\n";
}

// —— тестовая функция  ——————————————————————————

TestRunResult run_single_test_config(const TestConfig& cfg,
    int methodID /* 0=MDecomp,1=Rec,2=Simple */)
{
    const std::string dir = "C:/Users/User/source/repos/graduate work/graduate work/";
    std::vector<double> times;
    std::vector<double> reliabs;
    std::vector<long long> recs;

    for (int i = 0; i < cfg.num_repetitions; ++i) {
        std::ifstream fin(dir + cfg.graph_file_name);
        kGraph G = kGraphFileInput(fin);

        //G.CutPointsSearch(52, -1); //в переменной cutpoints будет лежать список точек сочленения

        // Сброс состояния, если нужно
        fin.close();

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

    const std::string dir = "C:/Users/User/source/repos/graduate work/graduate work/";

    std::vector<TestConfig> suite;
    for (const auto& pair : fileDiameters) {
        const std::string& fname = pair.first;
        const std::vector<int>& diamVec = pair.second;

        int s = -1, t = -1;
        {
            std::ifstream fin(dir + fname);
            if (!fin) {
                std::cerr << "Warning: cannot open " << fname << "\n";
                continue;
            }
            std::string line;
            std::getline(fin, line);
            std::getline(fin, line);
            if (std::getline(fin, line)) {
                std::istringstream ss(line);
                int val, idx = 0;
                while (ss >> val) {
                    if (val == 1) {
                        if (s < 0) s = idx;
                        else if (t < 0) t = idx;
                    }
                    if (ss.peek() == ',') ss.get();
                    ++idx;
                }
            }
            fin.close();
        }

        if (s < 1 || t < 1) {
            std::cerr << "Warning: no valid s/t found in " << fname << "\n";
            continue;
        }

        for (int D : diamVec) {
            suite.push_back({ fname, s, t, D, NUMBER_OF_REPETITIONS });
        }
    }

    std::ofstream sumF(dir + "test_summary_results.csv");
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