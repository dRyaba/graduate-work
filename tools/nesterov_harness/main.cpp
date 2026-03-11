/**
 * Nesterov CPFM harness - compare with graduate-work Method 4
 * K4, K5, 2_3x3_blocks, 3_blocks_sausage_3x3, 3_blocks_sausage_4x4
 * Note: Nesterov uses p.size()<=D (max path length D-1 edges). We pass D+1 to match graduate-work diameter D.
 */
#include "NesterovGraph.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

void runTest(const char* name, int n, int m, const std::vector<int>& kao,
             const std::vector<int>& fo_1based, const std::vector<int>& terminals,
             int D, double R) {
    std::cout << "\n=== " << name << " (n=" << n << ", m=" << m << ", D=" << D << ", R=" << R << ") ===\n";
    std::vector<int> fo = fo_1based;
    NesterovGraph g(kao.data(), fo.data(), terminals.data(), n, m, D, R);
    g.facto(true, true, true);
    std::cout << "RL = " << std::fixed << std::setprecision(15) << g.getRL() << "\n";
    std::cout << "RU = " << g.getRU() << "\n";
    std::cout << "Recursions = " << g.getRecursionsCount() << "\n";
}

bool loadKAOFile(const std::string& path, std::vector<int>& kao, std::vector<int>& fo,
                std::vector<int>& terminals, double& R) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;
    size_t comment = line.find("//");
    if (comment != std::string::npos) line = line.substr(0, comment);
    kao.clear();
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) kao.push_back(std::stoi(tok));
    }
    if (!std::getline(f, line)) return false;
    fo.clear();
    ss.str(line); ss.clear();
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) fo.push_back(std::stoi(tok));
    }
    if (!std::getline(f, line)) return false;
    terminals.clear();
    ss.str(line); ss.clear();
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) terminals.push_back(std::stoi(tok));
    }
    int n = static_cast<int>(kao.size()) - 1;
    if (n > 0 && terminals.size() == (size_t)(n + 1)) {
        terminals.erase(terminals.begin());
    }
    if (terminals.size() != (size_t)n) return false;
    if (!std::getline(f, line)) return false;
    comment = line.find("//");
    if (comment != std::string::npos) line = line.substr(0, comment);
    size_t c = line.find(',');
    if (c != std::string::npos) line[c] = '.';
    R = std::stod(line);
    return true;
}

int main() {
    std::cout << "Nesterov CPFM Harness (compare with graduate-work Method 4)\n";

    std::vector<int> k4_kao = {0, 3, 6, 9, 12};
    std::vector<int> k4_fo = {2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3};
    std::vector<int> k4_terminals = {1, 0, 0, 1};
    runTest("K4 (s=0, t=3, D=3)", 4, 6, k4_kao, k4_fo, k4_terminals, 4, 0.9);

    std::vector<int> k5_kao = {0, 4, 8, 12, 16, 20};
    std::vector<int> k5_fo = {2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4};
    std::vector<int> k5_terminals = {1, 0, 0, 0, 1};
    runTest("K5 (s=0, t=4, D=4)", 5, 10, k5_kao, k5_fo, k5_terminals, 5, 0.9);

    std::vector<int> kao, fo, terms;
    double R;
    const char* base = "graphs_data/";
    { std::ifstream p(std::string(base) + "2_3x3_blocks_kao.txt"); if (!p) base = "../graphs_data/"; }
    if (loadKAOFile(std::string(base) + "2_3x3_blocks_kao.txt", kao, fo, terms, R)) {
        int n = static_cast<int>(kao.size()) - 1;
        int m = static_cast<int>(fo.size()) / 2;
        runTest("2_3x3_blocks (s=0, t=16, D=10)", n, m, kao, fo, terms, 11, R);
    } else {
        std::cerr << "Could not load 2_3x3_blocks_kao.txt\n";
    }
    if (loadKAOFile(std::string(base) + "3_blocks_sausage_3x3_kao.txt", kao, fo, terms, R)) {
        int n = static_cast<int>(kao.size()) - 1;
        int m = static_cast<int>(fo.size()) / 2;
        runTest("3_blocks_sausage_3x3 (s=0, t=28, D=12)", n, m, kao, fo, terms, 13, R);
    } else {
        std::cerr << "Could not load 3_blocks_sausage_3x3_kao.txt\n";
    }
    std::cout << "\n--- 3_blocks_sausage_4x4 skipped (run manually if needed) ---\n";

    std::cout << "\nDone. Note: Nesterov counts both contract+delete as recursions.\n";
    return 0;
}
