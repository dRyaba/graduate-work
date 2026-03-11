/**
 * Portable Nesterov myGraph - CPFM for comparison with graduate-work
 * Replaced GetTickCount with std::chrono
 */
#include "NesterovGraph.h"
#include <cstring>

NesterovGraph::NesterovGraph(const int* kao, const int* fo, const int* terminals, int n, int m, int D, double R)
    : R(R), D(D), ppInfo(nullptr), ppathsWhichIncludeLinks(nullptr), pathsCount(0),
      sortedPathsWhichIncludeLinks(nullptr), probabilityOfLinks(nullptr) {
    KAO = std::vector<int>(kao, kao + (n + 1));
    FO = std::vector<int>(fo, fo + (2 * m));

    for (auto i = FO.begin(); i != FO.end(); ++i)
        *i = *i - 1;

    for (int i = 0; i < n; ++i) {
        if (terminals[i] == 1)
            K.push_back(i);
    }
    k = static_cast<int>(K.size());
    outFile.open("output_nesterov.txt", std::ios::out);
}

NesterovGraph::~NesterovGraph() {
    if (ppInfo != nullptr) {
        for (size_t i = 0; i < pData.size(); ++i)
            delete[] ppInfo[i];
        delete[] ppInfo;
    }
    if (ppathsWhichIncludeLinks != nullptr) {
        for (int i = 0; i < edgesCount; ++i)
            delete[] ppathsWhichIncludeLinks[i];
        delete[] ppathsWhichIncludeLinks;
    }
    if (sortedPathsWhichIncludeLinks != nullptr)
        delete[] sortedPathsWhichIncludeLinks;
    if (probabilityOfLinks != nullptr)
        delete[] probabilityOfLinks;
    outFile.close();
}

void NesterovGraph::evulatePst() {
    nesterov_path path(1);
    for (int i = 0; i < k - 1; ++i) {
        for (int j = i + 1; j < k; ++j) {
            path[0] = K[i];
            link st(K[i], K[j]);
            pData.push_back(pathsBetween2Terminals(st));
            evulatePst(st, path);
        }
    }
}

void NesterovGraph::evulatePst(link st, nesterov_path p) {
    int last = p.back(), s = st.first, t = st.second;
    if (last == t) {
        pData.back().paths.push_back(p);
        pathsCount++;
        return;
    }
    if (p.size() <= (unsigned int)D) {
        p.push_back(0);
        for (int i = KAO[last]; i < KAO[last + 1]; ++i) {
            bool visited = false;
            for (size_t j = 0; j < p.size() - 1; ++j) {
                if (p[j] == FO[i]) { visited = true; break; }
            }
            if (!visited) {
                p.back() = FO[i];
                evulatePst(st, p);
            }
        }
    }
}

void NesterovGraph::evulatePInfo() {
    int n = static_cast<int>(pData.size());
    ppInfo = new int*[n];
    for (int i = 0; i < n; ++i) {
        int m = static_cast<int>(pData[i].paths.size());
        ppInfo[i] = new int[2 + 2 * m];
        ppInfo[i][0] = 0;
        ppInfo[i][1] = m;
        for (int j = 0; j < m; ++j) {
            int l = static_cast<int>(pData[i].paths[j].size());
            ppInfo[i][2 + 2 * j] = 1;
            ppInfo[i][2 + 2 * j + 1] = l - 1;
        }
    }
}

void NesterovGraph::evulatePathsST() {
    int allEdgesCount = static_cast<int>(FO.size()) / 2;
    ppathsWhichIncludeLinks = new int*[allEdgesCount];
    edgesCount = 0;
    int n = static_cast<int>(KAO.size());
    for (int i = 0; i < n - 2; ++i) {
        for (int j = KAO[i]; j < KAO[i + 1]; ++j) {
            int vertex = FO[j];
            if (i < vertex) {
                ppathsWhichIncludeLinks[edgesCount++] = new int[pathsCount * 2 + 1 + 2];
                int linkPathsCount = 0;
                addToPathsST(i, vertex, linkPathsCount, edgesCount);
                ppathsWhichIncludeLinks[edgesCount - 1][0] = linkPathsCount;
                ppathsWhichIncludeLinks[edgesCount - 1][1] = i;
                ppathsWhichIncludeLinks[edgesCount - 1][2] = vertex;
                if (linkPathsCount == 0) {
                    edgesCount--;
                    delete[] ppathsWhichIncludeLinks[edgesCount];
                    continue;
                }
                if (linkPathsCount != pathsCount * 2) {
                    int* buf = new int[linkPathsCount * 2 + 1 + 2];
                    std::memcpy(buf, ppathsWhichIncludeLinks[edgesCount - 1],
                               (linkPathsCount * 2 + 1 + 2) * sizeof(int));
                    delete[] ppathsWhichIncludeLinks[edgesCount - 1];
                    ppathsWhichIncludeLinks[edgesCount - 1] = buf;
                }
            }
        }
    }
}

void NesterovGraph::addToPathsST(int vertexFrom, int vertexTo, int& linkPathsCount, int&) {
    for (size_t i = 0; i < pData.size(); ++i) {
        for (size_t j = 0; j < pData[i].paths.size(); ++j) {
            nesterov_path& p = pData[i].paths[j];
            for (size_t kk = 0; kk < p.size() - 1; ++kk) {
                if ((p[kk] == vertexFrom && p[kk + 1] == vertexTo) ||
                    (p[kk] == vertexTo && p[kk + 1] == vertexFrom)) {
                    if (p.size() == 2 && linkPathsCount != 0) {
                        ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * linkPathsCount] =
                            ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * 0];
                        ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * linkPathsCount + 1] =
                            ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * 0 + 1];
                        ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * 0] = static_cast<int>(i);
                        ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * 0 + 1] = static_cast<int>(j);
                        linkPathsCount++;
                        continue;
                    }
                    ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * linkPathsCount] = static_cast<int>(i);
                    ppathsWhichIncludeLinks[edgesCount - 1][1 + 2 + 2 * linkPathsCount + 1] = static_cast<int>(j);
                    linkPathsCount++;
                }
            }
        }
    }
}

void NesterovGraph::sortPathsSt(bool oldReduceOn, bool newReduceOn, bool oldSortOn) {
    sortedPathsWhichIncludeLinks = new int[edgesCount];
    probabilityOfLinks = new double[edgesCount];
    int reducedLinksCount = 0;
    for (int i = 0; i < edgesCount; ++i) {
        sortedPathsWhichIncludeLinks[i] = i;
        probabilityOfLinks[i] = R;
    }
    if (oldReduceOn) {
        for (int i = 0; i < edgesCount; ++i) {
            if (ppathsWhichIncludeLinks[i][0] == -1) continue;
            int u = ppathsWhichIncludeLinks[i][1];
            int v = ppathsWhichIncludeLinks[i][2];
            reduce(u, i, reducedLinksCount);
            reduce(v, i, reducedLinksCount);
        }
    }
    if (newReduceOn)
        reducedLinksCount += newReduce(ppInfo);
    if (oldSortOn)
        sort_quick(sortedPathsWhichIncludeLinks, 0, edgesCount - 1);
}

int NesterovGraph::newReduce(int** pInfo) {
    int reducedLinksCount = 0;
    for (int i = 0; i < edgesCount - 1; ++i) {
        if (ppathsWhichIncludeLinks[i][0] != -1) {
            for (int j = i + 1; j < edgesCount; ++j) {
                int iSize = ppathsWhichIncludeLinks[i][0];
                int jSize = ppathsWhichIncludeLinks[j][0];
                if (iSize == jSize && std::memcmp(&ppathsWhichIncludeLinks[i][3], &ppathsWhichIncludeLinks[j][3],
                                                  iSize * 2 * sizeof(int)) == 0) {
                    int n = ppathsWhichIncludeLinks[j][0];
                    for (int k = 0; k < n; ++k) {
                        int st = ppathsWhichIncludeLinks[i][1 + 2 + 2 * k];
                        int pSt = ppathsWhichIncludeLinks[i][1 + 2 + 2 * k + 1];
                        pInfo[st][2 + 2 * pSt + 1]--;
                    }
                    ppathsWhichIncludeLinks[j][0] = -1;
                    probabilityOfLinks[i] *= probabilityOfLinks[j];
                    reducedLinksCount++;
                }
            }
        }
    }
    return reducedLinksCount;
}

void NesterovGraph::sort_quick(int* a, int l, int r) {
    int i = l, j = r, x = a[(l + r) / 2];
    do {
        while (ppathsWhichIncludeLinks[a[i]][0] > ppathsWhichIncludeLinks[x][0]) i++;
        while (ppathsWhichIncludeLinks[x][0] > ppathsWhichIncludeLinks[a[j]][0]) j--;
        if (i <= j) {
            std::swap(a[i], a[j]);
            i++; j--;
        }
    } while (i < j);
    if (l < j) sort_quick(a, l, j);
    if (i < r) sort_quick(a, i, r);
}

void NesterovGraph::reduce(int& u, int& i, int& reducedLinksCount) {
    while (needReduce(u)) {
        for (int j = 0; j < edgesCount; ++j) {
            if (ppathsWhichIncludeLinks[j][0] == -1 || j == i) continue;
            if (ppathsWhichIncludeLinks[j][1] == u)
                u = ppathsWhichIncludeLinks[j][2];
            else if (ppathsWhichIncludeLinks[j][2] == u)
                u = ppathsWhichIncludeLinks[j][1];
            else continue;
            int n = ppathsWhichIncludeLinks[j][0];
            for (int k = 0; k < n; ++k) {
                int st = ppathsWhichIncludeLinks[i][1 + 2 + 2 * k];
                int pSt = ppathsWhichIncludeLinks[i][1 + 2 + 2 * k + 1];
                ppInfo[st][2 + 2 * pSt + 1]--;
            }
            ppathsWhichIncludeLinks[j][0] = -1;
            probabilityOfLinks[i] *= probabilityOfLinks[j];
            reducedLinksCount++;
            break;
        }
    }
}

bool NesterovGraph::needReduce(int& u) {
    for (int x : K)
        if (x == u) return false;
    return (KAO[u + 1] - KAO[u]) == 2;
}

int** NesterovGraph::copyPInfo(int** pInfo) {
    int n = static_cast<int>(pData.size());
    int** result = new int*[n];
    for (int i = 0; i < n; ++i) {
        if (pInfo[i][0] == 1) {
            result[i] = new int[1];
            result[i][0] = 1;
        } else {
            int m = static_cast<int>(pData[i].paths.size());
            result[i] = new int[2 + 2 * m];
            std::memcpy(result[i], pInfo[i], (2 + 2 * m) * sizeof(int));
        }
    }
    return result;
}

void NesterovGraph::facto(bool oldReduceOn, bool newReduceOn, bool oldSortOn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    evulatePst();
    evulatePInfo();
    auto t1 = std::chrono::high_resolution_clock::now();
    evulatePathsST();
    sortPathsSt(oldReduceOn, newReduceOn, oldSortOn);
    RL = 0;
    RU = 1;
    listsCount = 0;
    recursionsCount = 0;
    facto(copyPInfo(ppInfo), 0, 0, 1);
    auto t2 = std::chrono::high_resolution_clock::now();
    (void)t1;
    (void)t2;
}

void NesterovGraph::facto(int** pInfo, int conPairs, int eIndex, long double Pl) {
    int** buf1PInfo = copyPInfo(pInfo);
    contractLink(pInfo, conPairs, eIndex, Pl);
    deleteLink(buf1PInfo, conPairs, eIndex, Pl);
    for (size_t i = 0; i < pData.size(); ++i)
        delete[] buf1PInfo[i];
    delete[] buf1PInfo;
}

void NesterovGraph::contractLink(int** pInfo, int conPairs, int eIndex, long double Pl) {
    recursionsCount++;
    if (eIndex == edgesCount) return;
    int n = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][0];
    while (n == -1 && eIndex < edgesCount - 1)
        n = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[++eIndex]][0];
    Pl *= probabilityOfLinks[sortedPathsWhichIncludeLinks[eIndex]];
    bool needBreak = false;
    for (int i = 0; i < n; ++i) {
        int st = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][1 + 2 + 2 * i];
        int pSt = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][1 + 2 + 2 * i + 1];
        if (pInfo[st][0] == 0 && pInfo[st][2 + 2 * pSt] == 1) {
            pInfo[st][2 + 2 * pSt + 1]--;
            if (pInfo[st][2 + 2 * pSt + 1] == 0) {
                conPairs++;
                pInfo[st][0] = 1;
                if (conPairs == (k * (k - 1) / 2)) {
                    RL += Pl * 1.0;
                    RU -= Pl * 0.0;
                    needBreak = true;
                    break;
                }
            }
        }
    }
    if (!needBreak)
        facto(pInfo, conPairs, eIndex + 1, Pl);
}

void NesterovGraph::deleteLink(int** pInfo, int conPairs, int eIndex, long double Pl) {
    recursionsCount++;
    if (eIndex == edgesCount) return;
    int n = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][0];
    while (n == -1 && eIndex < edgesCount - 1)
        n = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[++eIndex]][0];
    Pl *= (1 - probabilityOfLinks[sortedPathsWhichIncludeLinks[eIndex]]);
    bool needBreak = false;
    for (int i = 0; i < n; ++i) {
        int st = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][1 + 2 + 2 * i];
        int pSt = ppathsWhichIncludeLinks[sortedPathsWhichIncludeLinks[eIndex]][1 + 2 + 2 * i + 1];
        if (pInfo[st][0] == 0 && pInfo[st][2 + 2 * pSt] == 1) {
            pInfo[st][2 + 2 * pSt] = 0;
            pInfo[st][1]--;
            if (pInfo[st][1] == 0) {
                RL += Pl * 0.0;
                RU -= Pl * 1.0;
                needBreak = true;
                break;
            }
        }
    }
    if (!needBreak)
        facto(pInfo, conPairs, eIndex + 1, Pl);
}

int NesterovGraph::getRecursionsCount() { return recursionsCount; }
long double NesterovGraph::getRL() { return RL; }
long double NesterovGraph::getRU() { return RU; }
