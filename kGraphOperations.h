#include <utility>

#ifndef GRADUATE_WORK_KGRAPHOPERATIONS_H
#define GRADUATE_WORK_KGRAPHOPERATIONS_H

struct Graph {
    std::vector<unsigned short> KAO, FO;
    std::vector<double> PArray;
    Graph() = default;
    Graph(std::vector<unsigned short> KAO, std::vector<unsigned short> FO, std::vector<double> PArray) {
        this->KAO = std::move(KAO);
        this->FO = std::move(FO);
        this->PArray = std::move(PArray);
    }
//    ~Graph()= default;
};

struct kGraph: public Graph {
    std::vector<unsigned short> Targets;

    kGraph() = default;

    kGraph(std::vector<unsigned short> KAO, std::vector<unsigned short> FO, std::vector<double> PArray,
           std::vector<unsigned short> Targets) : Graph() {
        this->KAO = std::move(KAO);
        this->FO = std::move(FO);
        this->PArray = std::move(PArray);
        this->Targets = std::move(Targets);
    }
//    ~kGraph()= default;
};



kGraph InducedKGraph(kGraph &G, std::vector<unsigned short> spisok, unsigned short Number);

int SearchEdge(const kGraph &G, int i, int j) {
//вычисляет номер ребра из i в j в массиве FO
    int res = G.FO.size();
    for (int k = G.KAO[i - 1]; k <= G.KAO[i] - 1; k++)
        if (G.FO[k] == j)
            res = k;
    return res;
}

kGraph DeleteEdgeK(const kGraph& G, unsigned short u, unsigned short v) {
//    Удаляет из графа ребро(u, v)
    unsigned short i, e, f;
    e = SearchEdge(G, u, v);
    f = SearchEdge(G, v, u);
    kGraph Result;
    Result.Targets.resize(G.Targets.size());
    if ((e < G.FO.size()) and (f < G.FO.size())) {
        Result.KAO.resize(G.KAO.size());
        Result.FO.resize(G.FO.size() - 2);
        Result.PArray.resize(G.PArray.size() - 2);
        if (e < f) {
            for (i = 0; i <= u - 1; i++)
                Result.KAO[i] = G.KAO[i];
            for (i = u; i <= v - 1; i++)
                Result.KAO[i] = G.KAO[i] - 1;
            for (i = v; i < G.KAO.size() - 1; i++)
                Result.KAO[i] = G.KAO[i] - 2;
            if (e > 0)
                for (i = 0; i <= e - 1; i++) {
                    Result.FO[i] = G.FO[i];
                    Result.PArray[i] = G.PArray[i];
                }
            if (f > e + 1)
                for (i = e; i < f - 2; i++) {
                    Result.FO[i] = G.FO[i + 1];
                    Result.PArray[i] = G.PArray[i + 1];
                }
            if (G.FO.size() > f + 1)
                for (i = f - 1; i < G.FO.size() - 3; i++) {
                    Result.FO[i] = G.FO[i + 2];
                    Result.PArray[i] = G.PArray[i + 2];
                }
        } else {
            for (i = 0; i <= v - 1; i++)
                Result.KAO[i] =G.KAO[i];
            for (i = v; i <= u - 1; i++)
                Result.KAO[i]=G.KAO[i] - 1;
            for (i = u; i <= G.KAO.size() - 1; i++)
                Result.KAO[i] = G.KAO[i] - 2;
            if (f > 0)
                for (i = 0; i <= f - 1; i++) {
                    Result.FO[i] = G.FO[i];
                    Result.PArray[i] = G.PArray[i];
                }
            if (e > f + 1)
                for (i = f; i <= e - 2; i++) {
                    Result.FO[i] = G.FO[i + 1];
                    Result.PArray[i] = G.PArray[i + 1];
                }
            if (G.FO.size() > e + 1)
                for (i = e - 1; i <= G.FO.size() - 3; i++) {
                    Result.FO[i] = G.FO[i + 2];
                    Result.PArray[i] = G.PArray[i + 2];
                }
        }
    }
    Result.Targets.resize(G.Targets.size());
    for (i = 0; i <= Result.Targets.size() - 1; i++)
        Result.Targets[i] = G.Targets[i];
    return Result;
}

bool KComponent(kGraph G) {
//{Выделяет компоненту связности графа G, содержащую все целевые вершины.
//Если это возможно, то компонента сохраняется как G, а результат - истина.
//В противном случае граф не изменяется, а результат - ложь.
//Переформирование графа G происходит лишь в том случае, когда он несвязен.
//Если граф не содержит целевых вершин, то результат - ложь.
//Преименение функции к пустому графу некорректно}
    int i, j, sum, l, sum1, SumAll;
    bool Result;
    std::vector<unsigned short> A, B, Spot;
    sum1 = 0;
    for (i = 1; i < G.Targets.size() - 1; i++)
        if (G.Targets[i] == 1)
            sum1++;
    if (sum1 == 1) {
        G.KAO.resize(2);
        G.KAO[0] = 0;
        G.KAO[1] = 0;
        G.FO.resize(0);
        G.PArray.resize(0);
        G.Targets.resize(2);
        G.Targets[0] = 0;
        G.Targets[1] = 0;
        Result = true;
    } else {
        sum = 1;
        SumAll = 1;
        A.resize(1);
        for (i = 1; i <= G.Targets.size() - 1; i++)
            if (G.Targets[i] == 1) {
                A[0] = i;
                break;
            }
        l = A.size();
        Spot.resize(G.KAO.size());
        for (i = 1; i < G.KAO.size() - 1; i++)
            Spot[i] = 2;
        Spot[A[0]] = 1;
        while (l > 0) {
            B.resize(0);
            for (i = 0; i < A.size() - 1; i++)
                for (j = G.KAO[A[i] - 1]; j < G.KAO[A[i] - 1]; j++)
                    if (Spot[G.FO[j]] == 2) {
                        B.resize(B.size() + 1);
                        B[B.size() - 1] = G.FO[j];
                        Spot[G.FO[j]] = 1;
                        if (G.Targets[G.FO[j]] == 1)
                            sum++;
                        SumAll++;
                    }

            A.resize(B.size());
            l = A.size();
            if (l > 0)
                for (i = 0; i <= l - 1; i++)
                    A[i] = B[i];
        }
        if (sum == sum1)
            Result = true;
        else
            Result = false;
        if (Result and (SumAll < (G.KAO.size() - 1)))
            G = InducedKGraph(G, Spot, 1);
    }
    return Result;
}

bool Boolka(kGraph &G, std::vector<unsigned short> spisok, unsigned short Number, unsigned short i, unsigned short j) {
    if (Number == 1)
        return true;
    else if ((spisok[i] == 0) and (spisok[G.FO[j]] == 0))
        return false;
    else
        return true;
}

kGraph InducedKGraph(kGraph &G, std::vector<unsigned short> spisok, unsigned short Number) {
//{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}
    std::vector<unsigned short> S, FO, KAO, Targets;
    unsigned short i, j, l;
    std::vector<double> PArray;
    S.resize(spisok.size());
    j = 0;
    for (i = 1; i <= spisok.size() - 1; i++)
        if ((spisok[i] = Number) or (spisok[i] = 0)) {
            j++;
            S[i] = j;
        } else
            S[i] = 0;
    KAO.resize(j + 1);
    KAO[0] = 0;
    FO.resize(G.FO.size());
    PArray.resize(G.FO.size());
    l = 0;
    for (i = 1; i <= spisok.size() - 1; i++)
        if (S[i] != 0) {
            KAO[S[i]] = KAO[S[i] - 1];
            for (j = G.KAO[i - 1]; j <= G.KAO[i - 1]; j++)
                if ((S[G.FO[j]] != 0) && Boolka(G, spisok, Number, i, j)) {
                    KAO[S[i]]++;
                    FO[l] = S[G.FO[j]];
                    PArray[l] = G.PArray[j];
                    l++;
                }
        }
    FO.resize(l);
    PArray.resize(l);
    Targets.resize(KAO.size());
    Targets[0] = 0;
    kGraph Result(KAO, FO, PArray, Targets);
    for (i = 1; G.Targets.size() - 1; i++)
        if (S[i] > 0)
            Result.Targets[S[i]] = G.Targets[i];

    return Result;
}
#endif //GRADUATE_WORK_KGRAPHOPERATIONS_H
