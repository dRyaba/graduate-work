/**
 * Portable copy of Nesterov myGraph for CPFM comparison.
 * Original: Nesterov_code/myGraphGUI/myGraph
 * Modified: removed Windows.h, use std::chrono for timing
 */
#ifndef NESTEROV_GRAPH_H
#define NESTEROV_GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <deque>
#include <chrono>

typedef std::vector<int> nesterov_path;

class NesterovGraph {
public:
    NesterovGraph(const int* kao, const int* fo, const int* terminals, int n, int m, int D, double R);
    ~NesterovGraph();

    void facto(bool oldReduceOn = false, bool newReduceOn = false, bool oldSortOn = false);
    int getRecursionsCount();
    long double getRL();
    long double getRU();

private:
    int k;
    int D;
    int recursionsCount;
    int listsCount;
    double R;

    std::vector<int> KAO;
    std::vector<int> FO;
    std::vector<int> K;

    typedef std::pair<int, int> link;
    struct pathsBetween2Terminals {
        link st;
        std::vector<nesterov_path> paths;
        pathsBetween2Terminals(link st) : st(st) {}
    };

    long double RL, RU;
    std::vector<pathsBetween2Terminals> pData;
    int** ppInfo;
    int** ppathsWhichIncludeLinks;
    int edgesCount;
    int pathsCount;
    int* sortedPathsWhichIncludeLinks;
    double* probabilityOfLinks;

    std::ofstream outFile;
    int min;

    int** copyPInfo(int** pInfo);
    void evulatePst();
    void evulatePst(link st, nesterov_path p);
    void evulatePInfo();
    void evulatePathsST();
    void addToPathsST(int vertexFrom, int vertexTo, int& linkPathsCount, int& edgesCount);
    void sortPathsSt(bool oldReduceOn, bool newReduceOn, bool oldSortOn);
    int newReduce(int** pInfo);
    void sort_quick(int* a, int l, int r);
    void reduce(int& u, int& i, int& reducedLinksCount);
    bool needReduce(int& u);
    void facto(int** pInfo, int conPairs, int eIndex, long double Pl);
    void contractLink(int** pInfo, int conPairs, int eIndex, long double Pl);
    void deleteLink(int** pInfo, int conPairs, int eIndex, long double Pl);
};

#endif
