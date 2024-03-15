#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <string>
#include "kGraphOperations.h"


void fileInput(std::vector<int>& KAO, std::vector<int>& FO, std::vector<int>& Targets, double& p) {
    std::ifstream fin("C:/Users/User/source/repos/graduate work/graduate work/input.txt");
    if (!fin) {
        std::cout << "Error!\n";
        throw std::runtime_error("OPEN_ERROR");
    }
    std::string line;
    std::getline(fin, line, '\n');
    std::istringstream skao(line);
    std::string temp;
    while (std::getline(skao, temp, ','))
        KAO.push_back(std::stoi(temp));

    std::getline(fin, line, '\n');
    std::istringstream sfo(line);
    while (std::getline(sfo, temp, ','))
        FO.push_back(std::stoi(temp));

    std::getline(fin, line, '\n');
    std::istringstream stargets(line);
    while (std::getline(stargets, temp, ','))
        Targets.push_back(std::stoi(temp));
    std::getline(fin, line);
    fin.close();
    size_t found = line.find(',');
    if (found != std::string::npos)
        line[found] = '.';
    p = std::stod(line);
}


int main() {
    std::cout << "Program is running\n";
    clock_t start_time = clock();
    clock_t total_start_time = start_time;
    std::vector<int> KAO, FO, Targets;
    double p;
    fileInput(KAO, FO, Targets, p);
    std::vector<double> Parray(FO.size(), p);
    kGraph G(KAO, FO, Parray, Targets);
    
    int x, y;
    for (int i = 0; i < Targets.size(); i++)
        if (Targets[i]) {
            x = i;
            break;
        }
    for (int i = x + 1; i < Targets.size(); i++)
        if (Targets[i]) {
            y = i;
            break;
        }
    int d = 8;
    G.ReliabilityDiamConstr2VertM(x, y, d);

    output << "Total Time(sec): " << (clock() - total_start_time) / 1000.0000;
    output.close();
    std::cout << "Program finished.";
}