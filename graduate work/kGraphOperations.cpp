#include <vector>
#include <fstream>
#include <iomanip>
#include "kGraphOperations.h"

std::ofstream output("C:/Users/User/source/repos/graduate work/graduate work/output.txt", std::ofstream::trunc);
int NumberOfRec = 0;
std::vector<double> sumReliab;
double globsumReliab = 0;


void Factoring(kGraph G, const int variant, const int d, double Reliab) {
    //результат лежит в sumReliab[0]
    //Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
    int i, j, k;
    NumberOfRec++;
    if (!variant) {
        if (!G.KComponent())
            return; // add 0 to sum
    }
    else if (!G.CheckDistanceFloyd(d))
        return; // add 0 to sum

    i = G.PArray.size();  //можно ли заменить PArray на FO? 
    for (j = G.PArray.size() - 1; j >= 0; j--)
        if (G.PArray[j] < 1) {
            i = j;
            break;
        }
    if (i == G.FO.size()) {
        sumReliab[0] += Reliab;
        return; // add Reliab to sum
    }

    double p = G.PArray[i];
    G.PArray[i] = 1;
    for (j = 1; j < G.KAO.size(); j++)
        if (G.KAO[j] > i)
            break;
    for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
        if (G.FO[k] == j)
            break;
    G.PArray[k] = 1;
    Factoring(G, 1, d, Reliab * p);

    kGraph T = G.DeleteEdgeK(G.FO[i], j);
    Factoring(T, 0, d, Reliab * (1 - p));
}

void Factoring2Vert(kGraph G, const int x, const int y, const int variant, const int d, double Reliab) {
    //результат лежит в sumReliab[0]
    //Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
    NumberOfRec++;
    if (!variant && G.DistanceDijkstra(x, y) > d)
        return;// add 0 to sum if distance from x to y >d
    int i = G.PArray.size();//можно ли заменить PArray на FO?
    for (int j = G.PArray.size() - 1; j >= 0; j--)
        if (G.PArray[j] < 1) {
            i = j;
            break;
        }
    if (i == G.FO.size()) {
        globsumReliab += Reliab;
        return; // add Reliab to sum
    }
    double p = G.PArray[i];
    G.PArray[i] = 1;
    int j;
    for (j = 1; j < G.KAO.size(); j++)
        if (G.KAO[j] > i)
            break;
    int k;
    for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
        if (G.FO[k] == j)
            break;
    G.PArray[k] = 1;
    Factoring2Vert(G, x, y, 1, d, Reliab * p);
    G.PArray[i] = p;
    G.PArray[k] = p;
    kGraph T = G.DeleteEdgeK(G.FO[i], j);
    Factoring2Vert(T, x, y, 0, d, Reliab * (1 - p));
}

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound) {
    //Заполняет вектор надежностей с соответствующим ограничением на диаметр
    //Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
    NumberOfRec++;
    int dist = G.DistanceDijkstra(x, y);
    if (!variant && dist > d) {
        if (dist <= UpperBound)
            Factoring2VertM(G, x, y, 0, d + 1, Reliab, LowerBound, UpperBound);
        return;// add 0 to sum if distance from x to y >d
    }
    int i = G.PArray.size();//можно ли заменить PArray на FO?
    for (int j = G.PArray.size() - 1; j >= 0; j--)
        if (G.PArray[j] < 1) {
            i = j;
            break;
        }
    if (i == G.FO.size()) {
        sumReliab[d - LowerBound] += Reliab;
        return; // add Reliab to sum
    }
    double p = G.PArray[i];
    G.PArray[i] = 1;
    int j;
    for (j = 1; j < G.KAO.size(); j++)
        if (G.KAO[j] > i)
            break;
    int k;
    for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
        if (G.FO[k] == j)
            break;
    G.PArray[k] = 1;
    Factoring2VertM(G, x, y, 1, d, Reliab * p, LowerBound, UpperBound);
    G.PArray[i] = p;
    G.PArray[k] = p;
    kGraph T = G.DeleteEdgeK(G.FO[i], j);
    Factoring2VertM(T, x, y, 0, d, Reliab * (1 - p), LowerBound, UpperBound);
}



void kGraph::ReliabilityDiamConstr(int d) {
    //{Расчет вероятности связности с ограничением на диаметр d методом ветвления.
    // Используется выделение комп-ты связности с целевыми в-ми, разделение ветвей,
    // встроенная ф-я проверки на расстояния}
    Nconst = this->KAO.size() * this->KAO.size();
    //double t = Factoring(G, 0, d);
    //return t;
    Factoring(*this, 0, d, 1);
}

void kGraph::ReliabilityDiamConstr2Vert(int x, int y, int d) {
    //    {Расчет в-ти св-ти двух вершин x,y с огр-м на диаметр d методом ветвления.
    //     Используется разделение ветвей, встроенная ф-я проверки на расстояния;
    //     не используется выделение компонентв с двумя в-ми (так быстрее)}
    Nconst = this->KAO.size() * this->KAO.size();
    Factoring2Vert(*this, x, y, 0, d, 1.0);
    output << std::setprecision(15) << globsumReliab << std::endl;
    output << "Recursions: " << NumberOfRec << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertM(int x, int y, int d) {
    //
    Nconst = this->KAO.size() * this->KAO.size();
    //std::vector<int> Blocks = G.DecomposeOnBlocks();

    int v = 9;
    int d1 = this->DistanceDijkstra(x, v), d2 = this->DistanceDijkstra(v, y);

    sumReliab.resize(d - d2 - d1 + 2);
    Factoring2VertM(*this, x, v, 0, d1 - 1, 1, d1 - 1, d - d2);
    std::vector<double> sumReliabG1(sumReliab);
    sumReliab.resize(0);
    sumReliab.resize(d - d1 - d2 + 1);
    Factoring2VertM(*this, v, y, 0, d2, 1, d2, d - d1);
    double totalRel = 0;
    for (int i = 0; i < d - d1 - d2 + 1; i++)
        totalRel += sumReliabG1[i + 1] * sumReliab[i];

    output << std::setprecision(15) << totalRel << std::endl;
    output << "Recursions: " << NumberOfRec << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertDecompose(int x, int y, int d) {
    // по идее формула используется та же, то есть нельзя ограничиваться одним значением d
    // для сетки 3х3 + 3х3, с объединением по одной вершине, если взять диаметр 8, то он будет единственным возможным
    //если взять диаметр 10, то там будет другой результат

    Nconst = this->KAO.size() * this->KAO.size();
    int v = 9;
    int d1 = this->DistanceDijkstra(x, v), d2 = this->DistanceDijkstra(v, y);

    Factoring2Vert(*this, x, v, 0, 4, 1);
    double totalRel = globsumReliab;
    globsumReliab = 0;
    Factoring2Vert(*this, v, y, 0, 4, 1);
    output << std::setprecision(15) << totalRel * globsumReliab << std::endl;
    output << "Recursions: " << NumberOfRec << std::endl;
}

bool kGraph::CheckDistanceFloyd(const int d) {
    //Методом Флойда строит матрицу расстояний и проверяет нужные
    int N = this->KAO.size();
    std::vector<std::vector<int> > M(N + 1);

    for (int i = 0; i < N; i++)
        M[i].resize(N + 1);
    for (int i = 1; i < N + 1; i++) {
        for (int j = i + 1; j < N + 1; j++) {
            if (this->CheckEdge(i, j))
                M[i][j] = 1;
            else M[i][j] = Nconst;
        }
        for (int k = 1; k < N + 1; k++) { //менял на +1
            for (i = 1; i < k; i++) {
                for (int j = i + 1; j < k; j++)
                    if (M[i][j] > M[i][k] + M[j][k])
                        M[i][j] = M[i][k] + M[j][k];
                for (int j = k + 1; j < N + 1; j++) //менял на +1
                    if (M[i][j] > M[i][k] + M[k][j])
                        M[i][j] = M[i][k] + M[k][j];
            }
            for (i = k + 1; i < N + 1; i++) { //менял на +1
                for (int j = i + 1; j < N + 1; j++) //менял на +1
                    if (M[i][j] > M[k][i] + M[k][j])
                        M[i][j] = M[k][i] + M[k][j];
            }
        }
    }
    for (int i = 1; i < N; i++) //убрал -1
        if (this->Targets[i])
            for (int j = i + 1; j < N; j++)
                if (this->Targets[j] && (M[i][j] > d))
                    return false;
    return true;
}

void kGraph::SearchVertice(std::vector<int>& DfNumber, std::vector<int>& TreeEdges, std::vector<int>& Low, std::vector<int>& Stek, int& r, int& l, int v, std::vector<int>& DOB) {
    int last = 0;
    DfNumber[v] = l;
    Low[v] = l;
    l++;
    for (int i = this->KAO[v - 1]; i < this->KAO[v]; i++) {
        if ((!DfNumber[this->FO[i]]) || ((DfNumber[this->FO[i]] < DfNumber[v]) && (!TreeEdges[i]))) {
            last = Stek.size();
            Stek.push_back(i);
            Stek.push_back(this->SearchEdge(this->FO[i], v));
        }
        if (!DfNumber[this->FO[i]]) {
            TreeEdges[i] = 1;
            TreeEdges[this->SearchEdge(this->FO[i], v)] = 1;
            this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, this->FO[i], DOB);

            if (Low[this->FO[i]] >= DfNumber[v]) {
                for (int j = last; j < Stek.size(); j++)
                    DOB[Stek[j]] = r;
                r++;
                Stek.resize(last);
            }
            Low[v] = std::min(Low[v], Low[this->FO[i]]);
            //if (Low[v] > Low[this->FO[i]]) 
            //	Low[v] = Low[this->FO[i]];
        }
        else if ((!TreeEdges[i]) && (DfNumber[this->FO[i]] < DfNumber[v]) && (Low[v] > DfNumber[this->FO[i]]))
            Low[v] = DfNumber[this->FO[i]];
    }
}

std::vector<int> kGraph::DecomposeOnBlocks() {
    //Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
    //Result[Result.size() - 1] - кол - во блоков
    std::vector<int> DfNumber(this->KAO.size()), TreeEdges(this->FO.size()), Low(this->KAO.size()), Stek(0);
    int r = 1, l = 1;
    std::vector<int> Result(this->FO.size() + 1);
    this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, 1, Result);
    Result[this->FO.size()] = r - 1;
    return Result;
}

kGraph kGraph::DeleteEdgeK(const int u, const int v) {
    //    Удаляет из графа ребро(u, v)
    int i, e, f;
    e = this->SearchEdge(u, v);
    f = this->SearchEdge(v, u);
    kGraph Result;
    Result.Targets.resize(this->Targets.size());
    if ((e < this->FO.size()) && (f < this->FO.size())) {
        Result.KAO.resize(this->KAO.size());
        Result.FO.resize(this->FO.size() - 2);
        Result.PArray.resize(this->PArray.size() - 2);
        if (e < f) {
            for (i = 0; i < u; i++)
                Result.KAO[i] = this->KAO[i];
            for (i = u; i < v; i++)
                Result.KAO[i] = this->KAO[i] - 1;
            for (i = v; i < this->KAO.size(); i++)
                Result.KAO[i] = this->KAO[i] - 2;
            if (e > 0)
                for (i = 0; i < e; i++) {
                    Result.FO[i] = this->FO[i];
                    Result.PArray[i] = this->PArray[i];
                }
            if (f > e + 1)
                for (i = e; i < f - 1; i++) {
                    Result.FO[i] = this->FO[i + 1];
                    Result.PArray[i] = this->PArray[i + 1];
                }
            if (this->FO.size() > f + 1)
                for (i = f - 1; i < this->FO.size() - 2; i++) {
                    Result.FO[i] = this->FO[i + 2];
                    Result.PArray[i] = this->PArray[i + 2];
                }
        }
        else {
            for (i = 0; i < v; i++)
                Result.KAO[i] = this->KAO[i];
            for (i = v; i < u; i++)
                Result.KAO[i] = this->KAO[i] - 1;
            for (i = u; i < this->KAO.size(); i++)
                Result.KAO[i] = this->KAO[i] - 2;
            if (f > 0)
                for (i = 0; i < f; i++) {
                    Result.FO[i] = this->FO[i];
                    Result.PArray[i] = this->PArray[i];
                }
            if (e > f + 1)
                for (i = f; i < e - 1; i++) {
                    Result.FO[i] = this->FO[i + 1];
                    Result.PArray[i] = this->PArray[i + 1];
                }
            if (this->FO.size() > e + 1)
                for (i = e - 1; i < this->FO.size() - 2; i++) {
                    Result.FO[i] = this->FO[i + 2];
                    Result.PArray[i] = this->PArray[i + 2];
                }
        }
    }
    Result.Targets.resize(this->Targets.size());
    for (i = 0; i < Result.Targets.size(); i++)
        Result.Targets[i] = this->Targets[i];
    return Result;
}

kGraph kGraph::InducedKGraph(std::vector<int> spisok, int Number) {
    //{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}
    std::vector<int> S, FO, KAO, Targets;
    std::vector<double> PArray;
    S.resize(spisok.size());
    int k = 0;
    for (int i = 1; i < spisok.size(); i++)
        if ((spisok[i] == Number) || (!spisok[i])) {
            k++;
            S[i] = k;
        }
        else
            S[i] = 0;
    KAO.resize(k + 1);
    KAO[0] = 0;
    FO.resize(this->FO.size());
    PArray.resize(this->FO.size());
    int l = 0;
    for (int i = 1; i < spisok.size(); i++)
        if (S[i] != 0) {
            KAO[S[i]] = KAO[S[i] - 1];
            for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
                if ((S[this->FO[j]]) && this->Boolka(spisok, Number, i, j)) {
                    KAO[S[i]]++;
                    FO[l] = S[this->FO[j]];
                    PArray[l] = this->PArray[j];
                    l++;
                }
        }
    FO.resize(l);
    PArray.resize(l);
    Targets.resize(KAO.size());
    Targets[0] = 0;
    kGraph Result(KAO, FO, PArray, Targets);
    for (int i = 1; i < this->Targets.size(); i++)
        if (S[i] > 0)
            Result.Targets[S[i]] = this->Targets[i];

    return Result;
}

bool kGraph::KComponent() {
    //{Выделяет компоненту связности графа G, содержащую все целевые вершины.
    //Если это возможно, то компонента сохраняется как G, а результат - истина.
    //В противном случае граф не изменяется, а результат - ложь.
    //Переформирование графа G происходит лишь в том случае, когда он несвязен.
    //Если граф не содержит целевых вершин, то результат - ложь.
    //Применение функции к пустому графу некорректно}
    int sum1 = 0;
    for (int i = 1; i < this->Targets.size(); i++)
        if (this->Targets[i] == 1)
            sum1++;
    if (sum1 == 1) {
        // как ищется компонента связности?
        // идем по массиву Таргетс со 2го до последнего элемента
        // и если там есть всего 1 целевая вершина, то возвращаем безвершинный граф
        // а если первый элемент таргетс 1
        this->KAO.resize(2);
        this->KAO[0] = 0;
        this->KAO[1] = 0;
        this->FO.resize(0);
        this->PArray.resize(0);
        this->Targets.resize(2);
        this->Targets[0] = 0;
        this->Targets[1] = 0;
        return true;
    }
    else {
        std::vector<int> A, B, Spot;
        int sum = 1;
        int SumAll = 1;
        A.resize(1);
        for (int i = 1; i < this->Targets.size(); i++)
            if (this->Targets[i] == 1) {
                A[0] = i;
                break;
            }
        int l = A.size();
        Spot.resize(this->KAO.size());
        for (int i = 1; i < this->KAO.size(); i++)
            Spot[i] = 2;
        Spot[A[0]] = 1;
        while (l > 0) {
            B.resize(0);
            for (int i : A)
                for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
                    if (Spot[this->FO[j]] == 2) {
                        B.resize(B.size() + 1);
                        B[B.size() - 1] = this->FO[j];
                        Spot[this->FO[j]] = 1;
                        if (this->Targets[this->FO[j]] == 1)
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
        if (Result && (SumAll < (this->KAO.size() - 1)))
            *this = this->InducedKGraph(Spot, 1); //requires attention!
        return Result;
    }
}

bool kGraph::Boolka(std::vector<int> spisok, int Number, int i, int j) {
    if (Number == 1)
        return true;
    //    if (!(spisok[i] || spisok[this->FO[j]]))
    //        return false;
    //    return true;
    return spisok[i] || spisok[this->FO[j]];
}