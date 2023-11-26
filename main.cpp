#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "kGraphOperations.h"


bool CheckEdge(kGraph &G, int i, int j);

bool CheckDistance(kGraph &G, int d, int Nconst);

double Factoring(kGraph &G, int variant, int NumberOfRec, int d, int Nconst);

double ReliabilityDiamConstr(kGraph &G, int d) {
//{Расчет вероятности связности с ограничением на диаметр d методом ветвления.
// Используется выделение комп-ты связности с целевыми в-ми, разделение ветвей,
// встроенная ф-я проверки на расстояния}
    int Nconst = G.KAO.size() * G.KAO.size();
    double t = Factoring(G, 0, 0, d, Nconst);
    return t;
}

double Factoring(kGraph &G, int variant, int NumberOfRec, int d, int Nconst) {
//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
    int i, j, k;
    NumberOfRec++;
    //правильно ли тут расставлены фигурные скобки
    if (variant == 0) {
        if (!KComponent(G))
            return 0;
    } else if (!CheckDistance(G, d, Nconst))
        return 0;

    i = G.PArray.size();
    for (j = G.PArray.size() - 1; j >= 0; j--)
        if (G.PArray[j] < 1) {
            i = j;
            break;
        }
    if (i == G.FO.size())
        return 1;
    else {
        double p = G.PArray[i];
        G.PArray[i] = 1;
        for (j = 1; j < G.KAO.size(); j++)
            if (G.KAO[j] > i)
                break;
        for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
            if (G.FO[k] == j)
                break;
        G.PArray[k] = 1;
        kGraph T = DeleteEdgeK(G, G.FO[i], j);
        double branch1 = p * Factoring(G, 1, NumberOfRec, d, Nconst);
        double branch2 = (1 - p) * Factoring(T, 0, NumberOfRec, d, Nconst);

        return branch1 + branch2;
    }
}

bool CheckDistance(kGraph &G, int d, int Nconst) {
//Методом Флойда строит матрицу расстояний и проверяет нужные
    int N = G.KAO.size();
    std::vector<std::vector<int> > M(N);

    for (int i = 0; i < N; i++) {
        M[i].resize(N + 1);
    }
    for (int i = 1; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (CheckEdge(G, i, j))
                M[i][j] = 1;
            else M[i][j] = Nconst;
        }
        for (int k = 1; k < N; k++) {
            for (i = 1; i < k; i++) {
                for (int j = i + 1; j < k; j++)
                    if (M[i][j] > M[i][k] + M[j][k])
                        M[i][j] = M[i][k] + M[j][k];
                for (int j = k + 1; j < N; j++)
                    if (M[i][j] > M[i][k] + M[k][j])
                        M[i][j] = M[i][k] + M[k][j];
            }
            for (i = k + 1; i < N; i++) {
                for (int j = i + 1; j < N; j++)
                    if (M[i][j] > M[k][i] + M[k][j])
                        M[i][j] = M[k][i] + M[k][j];
            }
        }
    }
    for (int i = 1; i < N - 1; i++)
        if (G.Targets[i] == 1)
            for (int j = i + 1; j < N - 1; j++)
                if ((G.Targets[j] == 1) && (M[i][j] > d))
                    return false;
    return true;
}

bool CheckEdge(kGraph &G, int i, int j) {
//вычисляет номер ребра из i в j в массиве FO
//ДР: скорее не вычисляет, а проверяет наличие
    for (int k = G.KAO[i - 1]; k < G.KAO[i]; k++)
        if (G.FO[k] == j)
            return true;
    return false;
}

int main() {
    std::ifstream fin("input.txt");
    if (!fin) {
        std::cout << "Error!" << std::endl;
        return 2;
    }
    std::vector<int> KAO, FO, Targets;

    std::string line;
    std::getline(fin, line,'\n');
    std::stringstream skao(line);
    std::string temp;
    while (getline(skao, temp, ','))
        KAO.push_back(std::stoi(temp));

    std::getline(fin, line,'\n');
    std::stringstream sfo(line);
    while (getline(sfo, temp, ','))
        FO.push_back(std::stoi(temp));

    std::getline(fin, line,'\n');
    std::stringstream stargets(line);
    while (getline(stargets, temp, ','))
        Targets.push_back(std::stoi(temp));
    std::getline(fin, line,'\n');
    fin.close();
    size_t found = line.find(',');
    if (found!=std::string::npos)
        line[found] = '.';
    double p = std::stod(line);
    int d = 4;
//    double p = 0.9;
//    KAO = {0, 4, 8, 12, 16, 20};
//    FO = {2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4};
//    Targets = {0, 1, 1, 1, 1, 0};
//    KAO = {0,1,3,4};
//    FO = {2,1,3,2};
//    Targets = {0,1,0,1};
    std::vector<double> Parray(FO.size(), p);
    kGraph G(KAO, FO, Parray, Targets);
    std::cout << std::setprecision(20) << ReliabilityDiamConstr(G, d);
}
