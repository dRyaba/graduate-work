#include <utility>

#ifndef GRADUATE_WORK_KGRAPHOPERATIONS_H
#define GRADUATE_WORK_KGRAPHOPERATIONS_H

struct Graph {
    std::vector<int> KAO, FO;
    std::vector<double> PArray;
    Graph() = default;
    Graph(std::vector<int> KAO, std::vector<int> FO, std::vector<double> PArray) {
        this->KAO = std::move(KAO);
        this->FO = std::move(FO);
        this->PArray = std::move(PArray);
    }
//    ~Graph()= default;
};

struct kGraph: public Graph {
    std::vector<int> Targets;

    kGraph() = default;

    kGraph(std::vector<int> KAO, std::vector<int> FO, std::vector<double> PArray,
           std::vector<int> Targets) : Graph() {
        this->KAO = std::move(KAO);
        this->FO = std::move(FO);
        this->PArray = std::move(PArray);
        this->Targets = std::move(Targets);
    }
//    ~kGraph()= default;
};



kGraph InducedKGraph(const kGraph &G, std::vector<int> spisok, int Number);

int SearchEdge(const kGraph &G, const int i,const int j) {
//вычисляет номер ребра из i в j в массиве FO
    int res = G.FO.size();
    for (int k = G.KAO[i - 1]; k < G.KAO[i]; k++)
        if (G.FO[k] == j)
            res = k;
    return res;
}

kGraph DeleteEdgeK(const kGraph& G,const int u,const int v) {
//    Удаляет из графа ребро(u, v)
    int i, e, f;
    e = SearchEdge(G, u, v);
    f = SearchEdge(G, v, u);
    kGraph Result;
    Result.Targets.resize(G.Targets.size());
    if ((e < G.FO.size()) && (f < G.FO.size())) {
        Result.KAO.resize(G.KAO.size());
        Result.FO.resize(G.FO.size() - 2);
        Result.PArray.resize(G.PArray.size() - 2);
        if (e < f) {
            for (i = 0; i < u; i++)
                Result.KAO[i] = G.KAO[i];
            for (i = u; i < v; i++)
                Result.KAO[i] = G.KAO[i] - 1;
            for (i = v; i < G.KAO.size(); i++)
                Result.KAO[i] = G.KAO[i] - 2;
            if (e > 0)
                for (i = 0; i < e; i++) {
                    Result.FO[i] = G.FO[i];
                    Result.PArray[i] = G.PArray[i];
                }
            if (f > e + 1)
                for (i = e; i < f - 1; i++) {
                    Result.FO[i] = G.FO[i + 1];
                    Result.PArray[i] = G.PArray[i + 1];
                }
            if (G.FO.size() > f + 1)
                for (i = f - 1; i < G.FO.size() - 2; i++) {
                    Result.FO[i] = G.FO[i + 2];
                    Result.PArray[i] = G.PArray[i + 2];
                }
        } else {
            for (i = 0; i < v; i++)
                Result.KAO[i] =G.KAO[i];
            for (i = v; i < u; i++)
                Result.KAO[i]=G.KAO[i] - 1;
            for (i = u; i < G.KAO.size(); i++)
                Result.KAO[i] = G.KAO[i] - 2;
            if (f > 0)
                for (i = 0; i < f; i++) {
                    Result.FO[i] = G.FO[i];
                    Result.PArray[i] = G.PArray[i];
                }
            if (e > f + 1)
                for (i = f; i < e - 1; i++) {
                    Result.FO[i] = G.FO[i + 1];
                    Result.PArray[i] = G.PArray[i + 1];
                }
            if (G.FO.size() > e + 1)
                for (i = e - 1; i < G.FO.size() - 2; i++) {
                    Result.FO[i] = G.FO[i + 2];
                    Result.PArray[i] = G.PArray[i + 2];
                }
        }
    }
    Result.Targets.resize(G.Targets.size());
    for (i = 0; i < Result.Targets.size(); i++)
        Result.Targets[i] = G.Targets[i];
    return Result;
}

bool KComponent(kGraph &G) {
//{Выделяет компоненту связности графа G, содержащую все целевые вершины.
//Если это возможно, то компонента сохраняется как G, а результат - истина.
//В противном случае граф не изменяется, а результат - ложь.
//Переформирование графа G происходит лишь в том случае, когда он несвязен.
//Если граф не содержит целевых вершин, то результат - ложь.
//Применение функции к пустому графу некорректно}
    int sum1 = 0;
    for (int i = 1; i < G.Targets.size(); i++)
        if (G.Targets[i] == 1)
            sum1++;
    if (sum1 == 1) {
        // как ищется компонента связности?
        // идем по массиву Таргетс со 2го до последнего элемента
        // и если там есть всего 1 целевая вершина, то возвращаем безвершинный граф
        // а если первый элемент таргетс 1
        G.KAO.resize(2);
        G.KAO[0] = 0;
        G.KAO[1] = 0;
        G.FO.resize(0);
        G.PArray.resize(0);
        G.Targets.resize(2);
        G.Targets[0] = 0;
        G.Targets[1] = 0;
        return true;
    } else {
        std::vector<int> A, B, Spot;
        int sum = 1;
        int SumAll = 1;
        A.resize(1);
        for (int i = 1; i < G.Targets.size(); i++)
            if (G.Targets[i] == 1) {
                A[0] = i;
                break;
            }
        int l = A.size();
        Spot.resize(G.KAO.size());
        for (int i = 1; i < G.KAO.size(); i++)
            Spot[i] = 2;
        Spot[A[0]] = 1;
        while (l > 0) {
            B.resize(0);
            for (int i : A)
                for (int j = G.KAO[i - 1]; j < G.KAO[i]; j++)
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
                for (int i = 0; i < l; i++)
                    A[i] = B[i];
        }
        bool Result = (sum == sum1);
        if (Result && (SumAll < (G.KAO.size() - 1)))
            G = InducedKGraph(G, Spot, 1);
        return Result;
    }
}

bool Boolka(const kGraph &G, std::vector<int> spisok, int Number, int i, int j) {
    if (Number == 1)
        return true;
//    if (!(spisok[i] || spisok[G.FO[j]]))
//        return false;
//    return true;
    return spisok[i] || spisok[G.FO[j]];
}

kGraph InducedKGraph(const kGraph &G, std::vector<int> spisok, int Number) {
//{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}
    std::vector<int> S, FO, KAO, Targets;
    std::vector<double> PArray;
    S.resize(spisok.size());
    int k = 0;
    for (int i = 1; i < spisok.size(); i++)
        if ((spisok[i] == Number) || (spisok[i] == 0)) {
            k++;
            S[i] = k;
        } else
            S[i] = 0;
    KAO.resize(k + 1);
    KAO[0] = 0;
    FO.resize(G.FO.size());
    PArray.resize(G.FO.size());
    int l = 0;
    for (int i = 1; i < spisok.size(); i++)
        if (S[i] != 0) {
            KAO[S[i]] = KAO[S[i] - 1];
            for (int j = G.KAO[i - 1]; j < G.KAO[i]; j++)
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
    for (int i = 1; i < G.Targets.size(); i++)
        if (S[i] > 0)
            Result.Targets[S[i]] = G.Targets[i];

    return Result;
}
#endif //GRADUATE_WORK_KGRAPHOPERATIONS_H
