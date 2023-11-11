//
// Created by User on 21.10.2023.
//

#ifndef GRADUATE_WORK_MATRIX_H
#define GRADUATE_WORK_MATRIX_H
#include <iostream>

class Matrix {
    int rows;
    int columns;
    int **values;
public:
    Matrix();

    Matrix(int rows, int columns);

    Matrix(int rows, int columns, int value);


    Matrix(int rows, int columns, int *values);

    Matrix(int rows, int columns, int **values);

// Получение минора
    Matrix getMinor(int a, int b);

// Сложение матриц
    Matrix operator+(const Matrix &a) const;

// Сравнение матриц
    bool operator==(const Matrix &b);

// Вывод матрицы на экран
    void operator+();
// Умножение матриц
    Matrix operator*(const Matrix &b) const ;

// Взятие строки
    int *operator[](int num);

// Взятие столбца
    int *operator()(int num);

// Транспонирование матрицы
    Matrix operator~() ;

    ~Matrix();
};
#endif //GRADUATE_WORK_MATRIX_H
