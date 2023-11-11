#include <iostream>
#include <vector>
#include "kGraphOperations.h"
#include "Matrix.h"

//проверить быстрее ли вычисление на статических массивах(мб ещё на calloc, malloc)

Matrix::Matrix() {
    this->rows = 0;
    this->columns = 0;
    this->values = nullptr;
}

Matrix::Matrix(int rows, int columns) {
    this->values = new int* [rows];
    for (int i = 0; i < rows; i++) {
        this->values[i] = new int[columns];
        for (int j = 0; j < columns; ++j)
            this->values[i][j] = 0;
        this->values[i][i] = 1;
    }
}

Matrix::Matrix(int rows, int columns, int value) : rows(rows), columns(columns) {
    this->values = new int* [rows];
    for (int i = 0; i < rows; i++) {
        this->values[i] = new int[columns];
        for (int j = 0; j < columns; ++j)
            this->values[i][j] = 0;
        this->values[i][i] = value;
    }
}

Matrix::Matrix(int rows, int columns, int *values) : rows(rows), columns(columns) {
    this->values = new int *[rows];
    for (int i = 0; i < rows; i++) {
        this->values[i] = new int[columns];
        for (int j = 0; j < columns; ++j)
            this->values[i][j] = 0;
        this->values[i][i] = values[i];
    }
}

Matrix::Matrix(int rows, int columns, int **values) : rows(rows), columns(columns) {
    this->values = values;
}

// Получение минора
Matrix Matrix::getMinor(int a, int b) {
    int **newValues = new int *[rows - 1];
    int ni = 0, nj = 0;
    for (int i = 0; i < rows; ++i) {
        if (i == a - 1)
            continue;
        newValues[ni] = new int[columns - 1];
        for (int j = 0; j < columns; ++j) {
            if (j == b - 1)
                continue;
            newValues[ni][nj++] = values[i][j];
        }
        ni++;
        nj = 0;
    }
    return {rows - 1, columns - 1, newValues};
}

// Сложение матриц
Matrix Matrix::operator+(const Matrix &a) const {
    if ((rows != a.rows) && (columns != a.columns)) {
        std::cout << "Ошибка: размеры матриц не совпадают.";
        exit(1);
    }
    int **newValues = new int *[rows];
    for (int i = 0; i < rows; i++) {
        newValues[i] = new int[columns];
        for (int j = 0; j < columns; j++)
            newValues[i][j] = this->values[i][j] + a.values[i][j];
    }
    return {rows, columns, newValues};
}

// Сравнение матриц
bool Matrix::operator==(const Matrix &b) {
    if ((rows != b.rows) && (columns != b.columns))
        return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            if (this->values[i][j] != b.values[i][j])
                return false;
    return true;
}

// Вывод матрицы на экран
void Matrix::operator+() {
    for (int i = 0; i < this->rows; ++i) {
        for (int j = 0; j < this->columns; ++j) {
            printf("%d ", this->values[i][j]);
        }
        printf("\n");
    }
}

// Умножение матриц
Matrix Matrix::operator*(const Matrix &b) const {
    if (columns != b.rows) {
        std::cout << "Ошибка: количество столбцов первой матрицы не совпадает с количеством строк второй матрицы.";
        exit(1);
    }
    int **newValues = new int *[rows];
    for (int i = 0; i < rows; i++) {
        newValues[i] = new int[b.columns];
        for (int j = 0; j < b.columns; j++) {
            int tempSum = 0;
            for (int k = 0; k < columns; ++k) {
                tempSum += this->values[i][k] * b.values[k][j];
            }
            newValues[i][j] = tempSum;
        }
    }
    return {rows, b.columns, newValues};
}

// Взятие строки
int *Matrix::operator[](int num) {
    return values[num];
}


// Транспонирование матрицы
Matrix Matrix::operator~() {
    int **newValues = new int *[rows];
    for (int i = 0; i < rows; ++i) {
        newValues[i] = new int[columns];
        for (int j = 0; j < columns; ++j)
            newValues[i][j] = this->values[j][i];
    }
    return {columns, rows, newValues};
}

Matrix::~Matrix() {
    for (int i = 0; i < columns; ++i)
        delete[] values[i];
    delete[] values;
}


bool CheckEdge(kGraph &G, int i, int j);

bool CheckDistance(kGraph &G, unsigned short d, unsigned short Nconst);

double Factoring(kGraph &G, unsigned short variant, int NumberOfRec, unsigned short d, unsigned short Nconst);

double ReliabilityDiamConstr(kGraph &G, unsigned short d) {
//{Расчет вероятности связности с ограничением на диаметр d методом ветвления.
// Используется выделение комп-ты связности с целевыми в-ми, разделение ветвей,
// встроенная ф-я проверки на расстояния}
    unsigned short Nconst;
    Nconst = (G.KAO.size()) * (G.KAO.size());
    double t = Factoring(G, 0, 0, d, Nconst);
    return t;
}

//Зачем NumberofRec? чтобы знать какой шаг ветвеления?
double Factoring(kGraph &G, unsigned short variant, int NumberOfRec, unsigned short d, unsigned short Nconst) {
//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
    double p;
    int i, j, k;
    int numberOfRec = NumberOfRec + 1;
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
        p = G.PArray[i];
        G.PArray[i] = 1;
        for (j = 1; j <= G.KAO.size() - 1; j++)
            if (G.KAO[j] > i)
                break;
//        int g = G.KAO[G.FO[i] - 1];
//        int g_top = G.KAO[G.FO[i]] - 1;
        for (k = G.KAO[G.FO[i] - 1]; k <= G.KAO[G.FO[i]] - 1; k++)
            if (G.FO[k] == j)
                break;
        G.PArray[k] = 1;
        kGraph T = DeleteEdgeK(G, G.FO[i], j);
        double a1 = p * Factoring(G, 1, numberOfRec, d, Nconst), a2 = (1 - p) * Factoring(T, 0, numberOfRec, d, Nconst);
        if (a1 + a2 != 1) {
            std::cout <<"a";
        }
        return a1 + a2;
    }
}

bool CheckDistance(kGraph &G, unsigned short d, unsigned short Nconst) {
//Методом Флойда строит матрицу расстояний и проверяет нужные
    unsigned short N = G.KAO.size() - 1;
    Matrix M(N + 1, N + 1);
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            if (CheckEdge(G, i, j))
                M[i][j] = 1;
            else M[i][j] = Nconst;
        }
        for (int k = 1; k < N; k++) {
            for (i = 1; i <= k - 1; i++) {
                for (int j = i + 1; j <= k - 1; j++)
                    if (M[i][j] > M[i][k] + M[j][k])
                        M[i][j] = M[i][k] + M[j][k];

                for (int j = k + 1; j <= N; j++)
                    if (M[i][j] > M[i][k] + M[k][j])
                        M[i][j] = M[i][k] + M[k][j];
            }
            for (i = k + 1; i <= N; i++) {
                for (int j = i + 1; j <= N; j++)
                    if (M[i][j] > M[k][i] + M[k][j])
                        M[i][j] = M[k][i] + M[k][j];
            }
        }
    }
    bool Result = true;
    for (int i = 1; i <= N - 1; i++)
        if (G.Targets[i] == 1)
            for (int j = i + 1; j < N; j++)
                if ((G.Targets[j] == 1) && (M[i][j] > d))
                    Result = false;
    return Result;
}

bool CheckEdge(kGraph& G, int i, int j) {
//вычисляет номер ребра из i в j в массиве FO
    bool Result = false;
    for (int k = G.KAO[i - 1]; k < G.KAO[i] - 1; k++)
        if (G.FO[k] == j)
            Result = true;
    return Result;
}

int main() {
    std::vector<unsigned short> KAO, FO, Targets;
    int d = 4;
    KAO = {0,1,3,4};
    FO = {2,1,3,2};
    std::vector<double> Parray(FO.size(), 0.9);
    Targets = {1,0,1};
//    Parray = {0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9};
    kGraph G(KAO,FO,Parray, Targets);
    std::cout<< ReliabilityDiamConstr(G,d);
}
