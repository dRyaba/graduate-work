#include <vector>
#include "GraphOperations.h"


int Graph::SearchEdge(const int i, const int j) {
    //вычисляет номер ребра из i в j в массиве FO
    int res = this->FO.size();
    for (int k = this->KAO[i - 1]; k < this->KAO[i]; k++)
        if (this->FO[k] == j)
            res = k;
    return res;
}

int Graph::DistanceDijkstra(const int x, const int y) {
    //    {Методом Дейкстры строит массив расстояний и возвращает нужное.
    //    Вначале - проверка вершины на изолированность}

    if (this->KAO[x] == this->KAO[x - 1]) //x не может быть 0 по условию задания KAO
        return Nconst;

    std::vector<int> DS(this->KAO.size(), Nconst), Pon(this->KAO.size());

    DS[x] = 0;
    int c, b = 0;
    for (int k = 0; k < this->KAO.size() - 2; k++) {
        c = Nconst;
        for (int i = 1; i < this->KAO.size(); i++)
            if ((DS[i] < c) && !Pon[i]) {
                c = DS[i];
                b = i;
            }
        Pon[b] = 1;
        for (int i = this->KAO[b - 1]; i < this->KAO[b]; i++) { // смотрим все ребра вершины b
            c = this->FO[i];
            if (!Pon[c] && (DS[c] > DS[b] + 1)) // ищем наименьший путь в вершину
                DS[c] = DS[b] + 1;
            //  {если длина ребра >1, то вместо +1 нужно +d(b,c)}
        }
    }
    return DS[y];
}

bool Graph::CheckEdge(int i, int j) {
    //вычисляет номер ребра из i в j в массиве FO
    //ДР: скорее не вычисляет, а проверяет наличие
    for (int k = this->KAO[i - 1]; k < this->KAO[i]; k++)
        if (this->FO[k] == j)
            return true;
    return false;
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
                    B[B.size() - 1] = this->FO[j];
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
    // {Находит минимальное сечение, если оно есть, то Result[0] = 2,
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
    //{Находит первое 2-сечение, если оно есть, то Result[0]=2, если нет - Result[0]=1;
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
                        Snew[Snew.size() - 1] = this->FO[j1];
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
