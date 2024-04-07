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
        return;// add 0 to sum if distance from x to y > d
    
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
        if (dist <= UpperBound){
            NumberOfRec--;
            Factoring2VertM(G, x, y, 0, d + 1, Reliab, LowerBound, UpperBound);
        }
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

kGraph UnionGraphs(kGraph G1, kGraph G2, int k) {
    //{Формирует граф, полученный объединением G1 и G2 по их первым k вершинам}
    //{NN - массив с номерами в-н G2 в G; номера в-н G1 в G совпадает с их номерами в G1}
    int N1 = G1.KAO.size() - 1, N2 = G2.KAO.size() - 1;
    std::vector<int> NN(N2 + 1);
    for (int i = 1; i < k + 1; i++)
        NN[i] = N1 + i - k;
    for (int i = k + 1; i < N2 + 1; i++)
        NN[i] = N1 + i - k;

    std::vector<int> KAO(N1 + N2 - k + 1), FO(0);
    int l = 0;
    bool boolean;
    for (int i = 1; i < k + 1; i++) {
        KAO[i] = KAO[i - 1];
        for (int j = G1.KAO[i - 1]; j < G1.KAO[i]; j++) {
            KAO[i]++;
            FO.resize(++l);
            FO[l - 1] = G1.FO[j];
        }
        for (int j = G2.KAO[i - 1]; j < G2.KAO[i]; j++) {
            boolean = true;
            for (int s = KAO[i - 1]; s < KAO[i]; s++)
                if (FO[s] == NN[G2.FO[j]])
                    boolean = false;
            if (boolean) {
                KAO[i]++;
                FO.resize(++l);
                FO[l - 1] = NN[G2.FO[j]];
            }
        }
    }
    for (int i = k + 1; i < N1 + 1; i++) {
        KAO[i] = KAO[i - 1];
        for (int j = G1.KAO[i - 1]; j < G1.KAO[i]; j++) {
            KAO[i]++;
            FO.resize(++l);
            FO[l - 1] = G1.FO[j];
        }
    }
    for (int i = N1 + 1; i < N1 + N2 - k + 1; i++) {
        KAO[i] = KAO[i - 1];
        for (int j = G2.KAO[i - N1 + k - 1]; j < G2.KAO[i - N1 + k]; j++) {
            KAO[i]++;
            FO.resize(++l);
            FO[l - 1] = NN[G2.FO[j]];
        }
    }
    std::vector<double> PArray(FO.size());
    for (int i = 0; PArray.size(); i++)
        PArray[i] = 0.9; // плохо, надо брать вероятности, соответствующие ребрам
    std::vector<int> Targets(KAO.size());
    //Targets тоже надо заполнять соответствующим образом
    kGraph Result(KAO, FO, PArray, Targets);

    return Result;
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

void kGraph::ReliabilityDiamConstr2VertM(int x, int y, const int& LowerBound, const int& UpperBound) {

    Nconst = this->KAO.size() * this->KAO.size();
    clock_t start_time = clock();
    this->CutPointsSearch(1, -1);
    if (cutPoints.empty()) {
        sumReliab.resize(UpperBound - LowerBound + 1);
        Factoring2VertM(*this, x, y, 0, UpperBound, 1, LowerBound, UpperBound);
        for (int i = 1; i < sumReliab.size(); i++)
            sumReliab[i] += sumReliab[i - 1];
        output << std::setprecision(15) << sumReliab[sumReliab.size() - 1] << std::endl;
        output << "Recursions: " << NumberOfRec << std::endl;
    }
    else {
        int v = cutPoints[0];
        int d1 = this->DistanceDijkstra(x, v), d2 = this->DistanceDijkstra(v, y);
        //int LowerBound = d1 - 1; //нужно ли тогда передавать эту переменную
        //int UpperBound = 10;
        sumReliab.resize(UpperBound - d1 - d2 + 2);

        //вручную создаём компоненту связности для графа 3х3
        kGraph T = *this;
        T.KAO.resize(v + 1);
        T.KAO[v] = 24;
        T.FO.resize(24);
        T.PArray.resize(24);
        T.Targets[v] = 1; T.Targets[T.Targets.size() - 1] = 0;

        //делаем для неё расчёт надёжности в заданном диаметре
        clock_t start_time = clock();
        Factoring2VertM(T, x, v, 0, d1, 1, d1 - 1, UpperBound - d2);
        output << "First part MFacto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
        output << "First part Recursions: " << NumberOfRec << std::endl;

        std::vector<double> sumReliabG1(sumReliab);
        ////костыль на получение надежности для отрезка диаметров
        //std::vector<double> ReliabDiam(sumReliab);
        //for (int i = 1; i < ReliabDiam.size(); i++)
        //    ReliabDiam[i] += ReliabDiam[i - 1];

        sumReliab.resize(0);
        sumReliab.resize(UpperBound - d1 - d2 + 1);

        //T = *this;
        //for (int i = 0; i < 24; i++){
        //    T.FO[i] = T.FO[i + 24] - 8;
        //    T.PArray[i] = T.PArray[i + 24];
        //}
        //T.FO.resize(24);
        //T.PArray.resize(24);
        //for (int i = 0; i < T.KAO.size(); i++) 
        //    T.KAO[i] = std::max(0, T.KAO[i] - 24);
        //for (int i = 1; i < 10; i++)
        //    T.KAO[i] = T.KAO[i + 8];
        //T.KAO.resize(10);


        start_time = clock();
        int NumOfRec1 = NumberOfRec;
        Factoring2VertM(T, x, v, 0, d2, 1, d2, UpperBound - d1);
        output << "Second part MFacto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
        output << "Second part Recursions: " << NumberOfRec - NumOfRec1 << std::endl;

        for (int i = 1; i < sumReliab.size(); i++)
            sumReliab[i] += sumReliab[i - 1];
        double totalRel = 0;
        for (int i = 1; i < UpperBound - d1 - d2 + 1; i++)
            totalRel += sumReliabG1[i] * sumReliab[UpperBound - d1 - d2 + 1 - i];

        output << std::setprecision(15) << totalRel << std::endl;
        output << "Recursions: " << NumberOfRec << std::endl;
    }
    output << "Modernized Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertDecompose(int x, int y, const int& LowerBound, const int& UpperBound) {
    Nconst = this->KAO.size() * this->KAO.size();
    clock_t start_time = clock();
    std::vector<int> spisok = this->DecomposeOnBlocksK();

    if (spisok[spisok.size() - 1] == 1) {
        sumReliab.resize(UpperBound - LowerBound + 1);
        Factoring2VertM(*this, x, y, 0, UpperBound, 1, LowerBound, UpperBound);
        for (int i = 1; i < sumReliab.size(); i++)
            sumReliab[i] += sumReliab[i - 1];
        output << std::setprecision(15) << sumReliab[sumReliab.size() - 1] << std::endl;
        output << "Recursions: " << NumberOfRec << std::endl;
        output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
        return;
    }

    kGraph Restored1 = this->RestoreBlock(1, spisok);
    for (int i = 0; i < Restored1.Targets.size(); i++)
        if (Restored1.Targets[i]) {
            x = i;
            break;
        }
    for (int i = x + 1; i < Restored1.Targets.size(); i++)
        if (Restored1.Targets[i]) {
            y = i;
            break;
        }
    int d1 = Restored1.DistanceDijkstra(x, y);
    kGraph Restored2 = this->RestoreBlock(2, spisok);
    //int LowerBound = d1 - 1;
    int x2, y2;
    for (int i = 0; i < Restored2.Targets.size(); i++)
        if (Restored2.Targets[i]) {
            x2 = i;
            break;
        }
    for (int i = x + 1; i < Restored2.Targets.size(); i++)
        if (Restored2.Targets[i]) {
            y2 = i;
            break;
        }
    int d2 = Restored2.DistanceDijkstra(x2, y2);
    sumReliab.resize(UpperBound - d1 - d2 + 2);

    /*clock_t */start_time = clock();
    Factoring2VertM(Restored1, x, y, 0, d1 - 1, 1, d1 - 1, UpperBound - d2 + 1);
    output << "First part MFacto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
    output << "First part Recursions: " << NumberOfRec << std::endl;

    std::vector<double> sumReliabG1(sumReliab);
    //костыль на получение надежности для отрезка диаметров
    std::vector<double> ReliabDiam(sumReliab);
    for (int i = 1; i < ReliabDiam.size(); i++)
        ReliabDiam[i] += ReliabDiam[i - 1];

        
        
        
    sumReliab.resize(0);
    sumReliab.resize(UpperBound - d1 - d2 + 1);

    start_time = clock();
    int NumOfRec1 = NumberOfRec;
    Factoring2VertM(Restored2, x2, y2, 0, d2, 1, d2, UpperBound - d1);
    output << "Second part MFacto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
    output << "Second part Recursions: " << NumberOfRec - NumOfRec1 << std::endl;

    for (int i = 1; i < sumReliab.size(); i++)
        sumReliab[i] += sumReliab[i - 1];
    double totalRel = 0;
    for (int i = d1; i < UpperBound - d2 + 1; i++)
    //for (int i = 1; i < UpperBound - d1 - d2 + 1; i++)
        totalRel += sumReliabG1[i - d1 + 1] * sumReliab[UpperBound - d1 - i];

    output << std::setprecision(15) << totalRel << std::endl;
    output << "Recursions: " << NumberOfRec << std::endl;
    output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
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
    Low[v] = l++;
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
    std::vector<int> DfNumber(this->KAO.size()), TreeEdges(this->FO.size()), Low(this->KAO.size()), Stek;
    int r = 1, l = 1;
    std::vector<int> Result(this->FO.size() + 1);
    this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, 1, Result);
    Result[this->FO.size()] = r - 1;
    return Result;
}

std::vector<int> kGraph::DecomposeOnBlocksK() {
    //{Разложение K - графа на блоки.
    //Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
    //Result[length(Result) - 1] - кол - во блоков, которые содержат цел.в - ны или являются связующими
    //Ребра, входящие в блоки, которые не содержат цел.в - н и не являются связующими входят блок с номером 0
    //Точки сочленения в связующих блоках заносятся в цлевые вершины графа G
    //Если в графе одна цел.в - на, то Result[i] = 0 для всех i
    //Если в графе нет ребер, то Result состоит из одного элемента, Result[0] = 0
    //Применение функции некорректно к несвязому графу и графу без целевых вершин}
    std::vector<int> S = this->DecomposeOnBlocks();
    if (S.size() == 1)
        return S;
    std::vector<int> TargetBlocks(S[S.size() - 1] + 1);
    for (int i = 1; i < this->Targets[this->Targets.size() - 1] + 1; i++)
        if (this->Targets[i] == 1) {
            bool boolean = true;
            int k = S[this->KAO[i - 1]];
            for (int j = this->KAO[i - 1] + 1; j < this->KAO[i]; j++)
                if (S[j] != k)
                    boolean = false;
            if (boolean)
                TargetBlocks[k] = 1;
        }

    for (int i = 1; i < TargetBlocks.size(); i++)
        if (!TargetBlocks[i] && !ConnectivityWithoutBlock(i, S))
            TargetBlocks[i] = 1;

    for (int i = 1; i < this->KAO.size(); i++)
        if (!this->Targets[i]) {
            bool boolean = true;
            int k = 0;
            for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
                if (TargetBlocks[S[j]] == 1)
                    if (!k)
                        k = S[j];
                    else if (S[j] != k)
                        boolean = false;
            if (!boolean)
                this->Targets[i] = 1;
        }
    std::vector<int> NewNumbers(TargetBlocks.size());
    for (int i = 0; i < NewNumbers.size(); i++)
        NewNumbers[i] = i;
    int k = 0;
    for (int i = 1; i < TargetBlocks.size(); i++)
        if (!TargetBlocks[i]) {
            k++;
            NewNumbers[i] = 0;
            for (int j = 1; j < NewNumbers.size(); j++)
                if (j > i)
                    NewNumbers[j] = NewNumbers[j] - 1;
        }
    for (int i = 0; i < S.size() - 1; i++) 
        S[i] = NewNumbers[S[i]];
    S[S.size() - 1] -= k;
    return S;
}

kGraph kGraph::DeleteEdgeK(const int u, const int v) {
    //    Удаляет из графа ребро(u, v)
    int e, f;
    e = this->SearchEdge(u, v);
    f = this->SearchEdge(v, u);
    kGraph Result;
    Result.Targets.resize(this->Targets.size());
    if ((e < this->FO.size()) && (f < this->FO.size())) {
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
    Result.Targets.resize(this->Targets.size());
    for (int i = 0; i < Result.Targets.size(); i++)
        Result.Targets[i] = this->Targets[i];
    return Result;
}

kGraph kGraph::InducedKGraph(std::vector<int> spisok, int Number) {
    //{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}
    std::vector<int> S(spisok.size());
    int k = 0;
    for (int i = 1; i < spisok.size(); i++)
        if ((spisok[i] == Number) || (!spisok[i])) {
            k++;
            S[i] = k;
        }
        else
            S[i] = 0;
    std::vector<int> KAO(k + 1), FO(this->FO.size());
    std::vector<double>PArray(this->FO.size());
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
    std::vector<int> Targets(KAO.size());
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

bool kGraph::ConnectivityWithoutBlock(int Block, std::vector<int> S) {
    std::vector<int> A(1);
    for (int i = 1; i < this->Targets.size(); i++)
        if (this->Targets[i] == 1) {
            A[0] = i;
            break;
        }
    int l = A.size(), sum = 1;
    std::vector<int> Spot(this->KAO.size()), B;
    Spot[A[0]] = 1;
    while (l > 0) {
        B.resize(0);
        for (int i = 0; i < A.size(); i++)
            for (int j = this->KAO[A[i] - 1]; j < this->KAO[A[i]]; j++)
                if ((!Spot[this->FO[j]]) && (S[j] != Block)) {
                    B.resize(B.size() + 1, this->FO[j]);
                    Spot[this->FO[j]] = 1;
                    sum += this->Targets[this->FO[j]] == 1;
                }
        A.resize(0);
        l = B.size();
        if (l > 0)
            for (int i = 0; i < l; i++)
                A.push_back(B[i]);
    }

    int sum1 = 0;
    for (int i = 1; i < this->Targets.size(); i++)
        sum += this->Targets[i] == 1;
    if (this->KAO.size() == 1)
        return true;
    return sum == sum1;
}

kGraph kGraph::RestoreBlockK(int Number, const std::vector<int>& spisok) {
    //{Восстанавливает блок K - графа G с номером Number,
    // используя spisok предварительно сформированный процедурой DecopmoseOnBlocksK
    // Эта же процедура уже пополнила спсиок цел.в - н необходимыми точками сочленения}
    std::vector<int> S(this->KAO.size()), BackS(this->KAO.size());
    int j = 0;
    for (int i = 0; i < spisok.size() - 1; i++)
        if ((spisok[i] == Number) && !S[this->FO[i]]) {
            S[this->FO[i]] = ++j;
            BackS[j] = this->FO[i];
        }
    std::vector<int> FO(this->FO.size()), KAO(j + 1), Targets(j + 1);
    std::vector<double> Parray(this->FO.size());
    BackS.resize(j + 1);
    for (int i = 1; i < Targets.size(); i++)
        Targets[i] = this->Targets[BackS[i]];
    int l = 0;
    for (int i = 1; i < BackS.size(); i++)
        for (j = this->KAO[BackS[i] - 1]; j < this->KAO[BackS[i]]; j++)
            if (spisok[j] == Number) {
                FO[l] = S[this->FO[j]];
                Parray[l] = this->PArray[j];
                KAO[i] = ++l;
            }
    FO.resize(l);
    Parray.resize(l);
    kGraph Result(KAO, FO, Parray, Targets);
    return Result;
}

kGraph kGraph::RestoreBlock(int Number, const std::vector<int>& spisok) {
    std::vector<int> S(this->KAO.size()), FO;
    std::vector<double> Parray(0);
    int j = 0; 
    int minFO = this->KAO.size();
    for (int i = 0; i < spisok.size() - 1; i++)
        if (spisok[i] == Number) {
            FO.push_back(this->FO[i]);
            Parray.push_back(this->PArray[i]);
            if (!S[this->FO[i]]){
                j++;
                minFO = std::min(minFO, this->FO[i]);
            }
            S[this->FO[i]]++;
            //BackS[j] = this->FO[i];
        }
    //std::vector<int> KAO(j + 1), Targets(j + 1);
    std::vector<int> KAO(1), Targets(1);


    for (int i = 1; i < S.size(); i++) {
        if (S[i]){
            KAO.push_back(KAO[KAO.size() - 1] + S[i]);
            Targets.push_back(this->Targets[i]);
        }
    }
    //for (int i = 1; i < KAO.size(); i++) {
    //    KAO[i] += S[i];
    //    if (S[i] && this->Targets[i])
    //        Targets[i] = 1;
    //}

    for (int i = 0; i < FO.size(); i++)
        FO[i] -= minFO - 1;
    kGraph Result(KAO, FO, Parray, Targets);
    return Result;
}