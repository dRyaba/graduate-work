#include <iostream>
#include <fstream>
//#include <sstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <string>
#include "kGraphOperations.h"




void GraphMerging(int k) {
    std::ifstream fin("C:/Users/User/source/repos/graduate work/graduate work/GraphsToMerge.txt");
    if (!fin) {
        std::cout << "Error!\n";
        throw std::runtime_error("OPEN_ERROR");
    }
    kGraph G1 = kGraphFileInput(fin);
    kGraph G2 = kGraphFileInput(fin);
    fin.close();
    kGraph G = UnionGraphs(G1, G2, k);
    for (int i = 1; i < G1.KAO.size() / 2; i++)
        G.ChangVertex(i, G1.KAO.size() - 1);
    std::ofstream fout("C:/Users/User/source/repos/graduate work/graduate work/MergedGraphs.txt", std::ofstream::trunc);
    if (!fout) {
        std::cout << "Error!\n";
        throw std::runtime_error("OPEN_ERROR");
    }
    G.kGraphFileOutput(fout);
}


int main() {
    std::cout << "Program is running\n";
    
    if (0) {
        GraphMerging(1);
        return 0;
    }
    clock_t total_start_time = clock();
    std::ifstream fin;
    fin.open("C:/Users/User/source/repos/graduate work/graduate work/input.txt", std::ios::in);
    kGraph G = kGraphFileInput(fin);
    fin.close();

    int x = 0, y = 0;
    for (int i = 0; i < G.Targets.size(); i++)
        if (G.Targets[i]) {
            x = i;
            break;
        }
    for (int i = x + 1; i < G.Targets.size(); i++)
        if (G.Targets[i]) {
            y = i;
            break;
        }
    int LowerBound = 4;
    int UpperBound = 12;
    
    //clock_t start_time = clock();
    //G.ReliabilityDiamConstr2Vert(x, y, UpperBound);
    //output << "Simple Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl << std::endl;
    //NumberOfRec = 0;
    //globsumReliab = 0;

    //G.ReliabilityDiamConstr2VertM(x, y, LowerBound, UpperBound);
    //sumReliab.resize(0);
    //NumberOfRec = 0;
    //output << std::endl;
    G.ReliabilityDiamConstr2VertDecompose(x, y, LowerBound, UpperBound);
    

    
    //sumReliab.resize(UpperBound - LowerBound + 1);

    ////делаем для неё расчёт надёжности в заданном диаметре
    //start_time = clock();
    //Factoring2VertM(G, x, y, 0, LowerBound, 1, LowerBound, UpperBound);
    //clock_t end_time = clock();


    ////получение надежности для отрезка диаметров
    //std::vector<double> ReliabDiam(sumReliab);
    //for (int i = 1; i < ReliabDiam.size(); i++)
    //    ReliabDiam[i] += ReliabDiam[i - 1];
    //output << std::setprecision(15) << ReliabDiam[ReliabDiam.size() - 1] << std::endl;
    //output << "Recursions: " << NumberOfRec << std::endl;
    //output << "MFacto for diameter Time(sec): " << (end_time - start_time) / 1000.0000 << std::endl;

    output << "Total Time(sec): " << (clock() - total_start_time) / 1000.0000;
    output.close();
    std::cout << "Program finished.";
}