#include <iostream>
#include <fstream>
//#include <sstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <string>
#include "kGraphOperations.h"


int main() {
    std::cout << "Program is running\n";

    //Раскомментировать при необходимости объединить графы. Исходные записываются в файле "GraphsToMerge", полученный в "MergedGraphs"
    //int k = 1;
    //GraphMerging(1);
    //std::cout << "Graphs are merged\n";
    //return 0;

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
    //задание верхней границы диаметра
    int UpperBound = 16;
    output << "Diameter: " << UpperBound << std::endl;
    //нет нумерации блоков в нужном порядке, 
    //поэтому для структур с количеством блоков >2 при неправильной изначальной нумерации блоков программа работать не будет(для GEANT'а они нумеруются вручную)

    //G.ReliabilityDiamConstr2Vert(x, y, UpperBound);
    //G.ReliabilityDiamConstr2VertDecomposeSimpleFacto(x, y, UpperBound);
    G.ReliabilityDiamConstr2VertMDecompose(x, y, UpperBound);
    //G.ReliabilityDiamConstr2VertMDecomposeParallel(x, y, UpperBound);

    //sumReliab.resize(UpperBound - LowerBound + 1);

    output << "Total Time(sec): " << (clock() - total_start_time) / 1000.0000;
    output.close();
    output1.close();
    std::cout << "Program finished.";
}