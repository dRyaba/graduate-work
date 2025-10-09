#include "GraphOperations.h"
#include <list>
#include <string>
#include <algorithm>
#include <sstream>
#include <set>
#include <queue>

int Nconst;

std::vector<int> cutPoints(0);

void Graph::CutPointsSearch(int v, int p) {
    static std::vector<bool> used(this->KAO.size(), false);
    static std::vector<int> tin(this->KAO.size()), fup(this->KAO.size());
    static int timer = 0;
    used[v] = true;
    tin[v] = fup[v] = timer++;
    int children = 0;
    bool cutpointFlag = false;

    for (size_t i = this->KAO[v - 1]; i < this->KAO[v]; ++i) {
        int to = this->FO[i];
        //if (to == p)  continue;
        if (used[to])
            fup[v] = std::min(fup[v], tin[to]);
        else {
            this->CutPointsSearch(to, v);
            fup[v] = std::min(fup[v], fup[to]);
            if (fup[to] >= tin[v] && p != -1)
                cutpointFlag = true;
            ++children;
        }
    }
    if ((p == -1 && children > 1) || cutpointFlag)
        cutPoints.push_back(v);
}

int Graph::SearchEdge(const int i, const int j) const {
    // ищем ровно одно ребро (i->j) в CSR-массиве FO
    int start = KAO[i - 1];
    int end = KAO[i];
    for (int idx = start; idx < end; ++idx) {
        if (FO[idx] == j)
            return idx;            // нашли Ч сразу возвращаем
    }
    return static_cast<int>(FO.size());  // не нашли Ч возвращаем Ђнедопустимыйї индекс
}

int Graph::DistanceDijkstra(int src, int dst) const {
    // ѕроверка изолированности: если у src нет соседей
    if (KAO[src] == KAO[src - 1])
        return Nconst;

    const int n = static_cast<int>(KAO.size());
    const int INF = Nconst;

    // ћин-куча пар (dist, vertex)
    using PII = std::pair<int, int>;
    std::priority_queue<PII, std::vector<PII>, std::greater<PII>> pq;

    std::vector<int> dist(n, INF);
    std::vector<bool> seen(n, false);

    dist[src] = 0;
    pq.emplace(0, src);

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (seen[u]) continue;
        seen[u] = true;
        if (u == dst)   // можем прервать, как только дошли до dst
            return d;

        // –ассматриваем всех соседей u
        for (int idx = KAO[u - 1]; idx < KAO[u]; ++idx) {
            int v = FO[idx];
            int w = 1; // тут у вас вес = 1; если не 1, замените на нужное
            if (dist[v] > d + w && !seen[v]) {
                dist[v] = d + w;
                pq.emplace(dist[v], v);
            }
        }
    }

    return dist[dst];
}


bool Graph::CheckEdge(int i, int j) const {
    //вычисл€ет номер ребра из i в j в массиве FO
    //ƒ–: скорее не вычисл€ет, а провер€ет наличие
    auto begin = FO.begin() + KAO[i-1];
    auto end   = FO.begin() + KAO[i];
    return std::find(begin, end, j) != end;
}

std::vector<int> Graph::CutDecompose(const std::vector<int> V) {
    std::vector<int> Spot(this->KAO.size());
    for (int i = 1; i < this->KAO.size(); i++)
        Spot[i] = 1;
    for (int i = 0; i < V.size(); i++)
        Spot[V[i]] = 0;
    int c;
    for (int i = 1; i < this->KAO.size(); i++)
        if (Spot[i] == 1) {
            c = i;
            break;
        }
    Spot[c] = 2;
    int sum = 1, l = 1;
    std::vector<int> A(1, c), B;
    while (l > 0) {
        B.resize(0);
        for (int i = 0; i < A.size(); i++)
            for (int j = this->KAO[A[i] - 1]; j < this->KAO[A[i]]; j++)
                if (Spot[this->FO[j]] == 1) {
                    B.resize(B.size() + 1);
                    B.back() = this->FO[j];
                    Spot[this->FO[j]] = 2;
                    sum++;
                }
        A.resize(B.size());
        l = A.size();
        if (l > 0)
            for (int i = 0; i < l; i++)
                A[i] = B[i];
    }
    if (sum + V.size() == this->KAO.size() - 1)
        Spot[0] = 1;
    else
        Spot[0] = 2;
    return Spot;
}

Graph Graph::DeleteEdge(int u, int v) {
    int e, f;
    e = this->SearchEdge(u, v);
    f = this->SearchEdge(v, u);
    Graph Result;
    if (e < this->FO.size() && f < this->FO.size()) {
        Result.KAO.resize(this->KAO.size());
        Result.FO.resize(this->FO.size() - 2);
        Result.PArray.resize(this->PArray.size() - 2);
        if (e < f) {
            for (int i = 0; i < u; i++)
                Result.KAO[i] = this->KAO[i];
            for (int i = u; i < v; i++)
                Result.KAO[i] = this->KAO[i] - 1;
            for (int i = v; i < this->KAO.size(); i++)
                Result.KAO[i] = this->KAO[i] - 2;
            if (e > 0)
                for (int i = 0; i < e; i++) {
                    Result.FO[i] = this->FO[i];
                    Result.PArray[i] = this->PArray[i];
                }
            if (f > e + 1)
                for (int i = e; i < f - 1; i++) {
                    Result.FO[i] = this->FO[i + 1];
                    Result.PArray[i] = this->PArray[i + 1];
                }
            if (this->FO.size() > f + 1)
                for (int i = f - 1; i < this->FO.size() - 2; i++) {
                    Result.FO[i] = this->FO[i + 2];
                    Result.PArray[i] = this->PArray[i + 2];
                }
        }
        else {
            for (int i = 0; i < v; i++)
                Result.KAO[i] = this->KAO[i];
            for (int i = v; i < u; i++)
                Result.KAO[i] = this->KAO[i] - 1;
            for (int i = u; i < this->KAO.size(); i++)
                Result.KAO[i] = this->KAO[i] - 2;
            if (f > 0)
                for (int i = 0; i < f; i++) {
                    Result.FO[i] = this->FO[i];
                    Result.PArray[i] = this->PArray[i];
                }
            if (e > f + 1)
                for (int i = f; i < e - 1; i++) {
                    Result.FO[i] = this->FO[i + 1];
                    Result.PArray[i] = this->PArray[i + 1];
                }
            if (this->FO.size() > e + 1)
                for (int i = e - 1; i < this->FO.size() - 2; i++) {
                    Result.FO[i] = this->FO[i + 2];
                    Result.PArray[i] = this->PArray[i + 2];
                }
        }
    }
    return Result;
}

std::vector<int> Graph::CutSearch() {
    // {Ќаходит минимальное сечение, если оно есть, то Result[0] = 2,
    // если нет - Result[0] = 1; Result[v] = 1, если v в G1, 
    // Result[v] = 2, если v в G2, Result[v] = 0, если v в сечении}
    std::vector<int> V, M;
    bool boolean;
    int N = this->KAO.size() - 1;
    for (int k = 2; k < N; k++) {
        V.resize(k);
        for (int i = 1; i < k + 1; i++)
            V[i - 1] = i;
        M = this->CutDecompose(V);
        boolean = (M[0] == 2);
        int l;
        while (boolean) {
            l = N + 1;
            for (int i = k - 1; i >= 0; i--)
                if (V[i] < N - (k - 1 - i)) {
                    l = i;
                    break;
                }
            if (l == N + 1)
                boolean = false;
            else {
                V[l]++;
                for (int i = l + 1; i < k; i++)
                    V[i] = V[i - 1] + 1;
                M = this->CutDecompose(V);
                if (M[0] == 2)
                    boolean = false;
            }
        }
        if (M[0] == 2)
            break;
    }
    return M;
}

std::vector<int> Graph::CutDecomposeOnTwo() {
    //{Ќаходит первое 2-сечение, если оно есть, то Result[0]=2, если нет - Result[0]=1;
    // Result[v] = 1, если v в G1, Result[v] = 2, если v в G2, Result[v] = 0, если v в сечении}
    bool boolean = true;
    int i = 1, j = 2, i1, j1, l, N;
    std::vector<int> Spot, S, Snew;
    N = this->KAO.size();
    while (boolean) {
        Spot.resize(N + 1);
        Spot[0] = 2;
        Spot[i] = 0;
        Spot[j] = 0;
        l = 3;
        S.resize(1);
        if (i > 1)
            i1 = 1;
        else if (j > 2)
            i1 = 2;
        else i1 = 3;
        S[0] = i1;
        Spot[i1] = 1;
        while (S.size() > 0) {
            Snew.resize(0);
            for (i1 = 0; i1 < S.size(); i1++)
                for (j1 = this->KAO[S[i1] - 1]; j1 < this->KAO[S[i1]]; j1++)
                    if (Spot[this->FO[j1]] == 2) {
                        Snew.resize(Snew.size() + 1);
                        Snew.back() = this->FO[j1];
                        Spot[this->FO[j1]] = 1;
                        l++;
                    }
            S.resize(Snew.size());
            if (S.size() > 0)
                for (i1 = 0; i1 < S.size(); i1++)
                    S[i1] = Snew[i1];
        }
        if (l == N)
            Spot[0] = 1;
        else
            boolean = false;
        if (j < N)
            j++;
        else if (i < N - 1) {
            i++;
            j++;
        }
        else
            boolean = false;
    }
    return Spot;
}

Graph Graph::ChangVertex(int u, int v) {
    //ћен€ет в графе вершины u и v местами (перенумеровывает)
    if (u > this->KAO.size() - 1 || v > this->KAO.size() - 1 || u == v)
        return *this;
    if (u > v)
        std::swap(u, v);
    std::vector<int> KAO(this->KAO.size()), FO(this->FO.size());
    std::vector<double> PArray(this->PArray.size());
    for (int i = 0; i < u; i++)
        KAO[i] = this->KAO[i];
    KAO[u] = this->KAO[u - 1] + this->KAO[v] - this->KAO[v - 1];
    for (int i = u + 1; i < v; i++)
        KAO[i] = KAO[i - 1] + this->KAO[i] - this->KAO[i - 1];
    KAO[v] = KAO[v - 1] + this->KAO[u] - this->KAO[u - 1];
    for (int i = v + 1; i < this->KAO.size(); i++)
        KAO[i] = this->KAO[i];

    int j = 0;
    for (int i = this->KAO[v - 1]; i < this->KAO[v]; i++) {
        FO[KAO[u - 1] + j] = this->FO[i];
        PArray[KAO[u - 1] + j] = this->PArray[i];
        j++;
    }
    j = 0;
    for (int i = this->KAO[u - 1]; i < this->KAO[u]; i++) {
        FO[KAO[v - 1] + j] = this->FO[i];
        PArray[KAO[v - 1] + j] = this->PArray[i];
        j++;
    }
    for (int ver = u + 1; ver < v; ver++) {
        j = 0;
        for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
            FO[KAO[ver - 1] + j] = this->FO[i];
            PArray[KAO[ver - 1] + j] = this->PArray[i];
            j++;
        }
    }
        
    for (int ver = 1; ver < u; ver++)
        for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
            FO[i] = this->FO[i];
            PArray[i] = this->PArray[i];
        }
    for (int ver = v + 1; ver < this->KAO.size(); ver++)
        for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
            FO[i] = this->FO[i];
            PArray[i] = this->PArray[i];
        }
    for (int i = 0; i < this->FO.size(); i++)
        if (FO[i] == u)
            FO[i] = v;
        else if (FO[i] == v)
            FO[i] = u;
    Graph Result(KAO, FO, PArray);
    return Result;
}

void Graph::convertEdgeListToKAOFO(const std::string& inputPath,
    const std::string& outputPath,
    double reliability) {
    std::ifstream inFile(inputPath);
    std::vector<std::pair<int, int>> edges;
    int u, v; char dash;
    while (inFile >> u >> dash >> dash >> v) {
        edges.emplace_back(u, v);
    }

    int maxVertex = 0;
    for (auto& e : edges) maxVertex = std::max(maxVertex, std::max(e.first, e.second));
    int n = maxVertex;

    std::vector<std::vector<int>> adj(n + 1);
    for (auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    for (int i = 1; i <= n; ++i) std::sort(adj[i].begin(), adj[i].end());

    KAO.resize(n + 2);
    KAO[1] = 0;
    for (int i = 1; i <= n; ++i) KAO[i + 1] = KAO[i] + static_cast<int>(adj[i].size());

    int totalSize = KAO[n + 1];
    FO.resize(totalSize);
    int idx = 0;
    for (int i = 1; i <= n; ++i)
        for (int neighbor : adj[i])
            FO[idx++] = neighbor;

    PArray.assign(totalSize, reliability);

    std::ofstream outFile(outputPath);
    for (int i = 1; i <= n + 1; ++i) outFile << KAO[i] << (i < n + 1 ? ',' : '\n');
    for (int i = 0; i < totalSize; ++i) outFile << FO[i] << (i + 1 < totalSize ? ',' : '\n');
    for (int i = 0; i < totalSize; ++i) outFile << PArray[i] << (i + 1 < totalSize ? ',' : '\n');
}

void Graph::convertKAOFOToEdgeList(const std::string& inputPath,
    const std::string& outputPath) {
    std::ifstream inFile(inputPath);
    std::string line, token;
    std::getline(inFile, line);
    std::stringstream ss(line);
    KAO.clear();
    while (std::getline(ss, token, ',')) KAO.push_back(std::stoi(token));

    std::getline(inFile, line);
    ss.clear(); ss.str(line);
    FO.clear();
    while (std::getline(ss, token, ',')) FO.push_back(std::stoi(token));

    std::ofstream outFile(outputPath);
    std::set<std::pair<int, int>> seen;
    int n = static_cast<int>(KAO.size()) - 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = KAO[i - 1]; j < KAO[i]; ++j) {
            int a = i, b = FO[j];
            if (a > b) std::swap(a, b);
            if (seen.insert({ a,b }).second)
                outFile << a << " -- " << b << '\n';
        }
    }
}
