//
// Created by Dmytro Klishch on 10/31/22.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <fstream>
#include <iostream>
#include "doctest.h"

void calculateLNonFirstColumn(int j, double **A, double **L, double **U, int n) {
    for (int i = j; i < n; ++i) {
        L[i][j] = A[i][j];
        for (int k = 0; k < j; ++k) {
            L[i][j] -= U[k][j] * L[i][k];
        }
    }
}

void calculateLFirstColumn(double **A, double **L, int n) {
    for (int i = 0; i < n; ++i) {
        L[i][0] = A[i][0];
    }
}

void calculateLColumn(int j, double **A, double **L, double **U, int n) {
    if (j == 0) calculateLFirstColumn(A, L, n);
    else calculateLNonFirstColumn(j, A, L, U, n);
}

void calculateUNonFirstRow(int i, double **A, double **L, double **U, int n) {
    for (int j = i; j < n; ++j) {
        U[i][j] = A[i][j];
        for (int k = 0; k < i; ++k) {
            U[i][j] -= L[i][k] * U[k][j];
        }
        U[i][j] /= L[i][i];
    }
}

void calculateUFirstRow(double **A, double **U, int n) {
    for (int j = 0; j < n; ++j) {
        U[0][j] = A[0][j] / A[0][0];
    }
}

void calculateURow(int i, double **A, double **L, double **U, int n) {
    if (i == 0) calculateUFirstRow(A, U, n);
    else calculateUNonFirstRow(i, A, L, U, n);
}

double **creteSquareMatrix(int n) {
    auto **matrix = new double *[n];
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = 0;
        }
    }
    return matrix;
}

void calculateLSubSolution(int k, double *d, double **L, int n) {
    d[0] = (k == 0) / L[0][0];
    for (int i = 1; i < n; ++i) {
        d[i] = (i == k);
        for (int j = 0; j < i; ++j) {
            d[i] -= L[i][j] * d[j];
        }
        d[i] /= L[i][i];
    }
}

void calculateSolution(int k, double **R, double *d, double **U, int n) {
    R[n - 1][k] = d[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        R[i][k] = d[i];
        for (int j = i + 1; j < n; ++j) {
            R[i][k] -= U[i][j] * R[j][k];
        }
    }
}

void calculateRColumn(int k, double **R, double **L, double **U, int n) {
    auto *d = new double[n];
    calculateLSubSolution(k, d, L, n);
    calculateSolution(k, R, d, U, n);
    delete[] d;
}

double **calculateInverseMatrix(double **A, int n) {
    double **L = creteSquareMatrix(n);
    double **U = creteSquareMatrix(n);
    double **R = creteSquareMatrix(n);

    for (int i = 0; i < n; ++i) {
        calculateLColumn(i, A, L, U, n);
        calculateURow(i, A, L, U, n);
    }

    for (int i = 0; i < n; ++i) {
        calculateRColumn(i, R, L, U, n);
    }

    delete[] L;
    delete[] U;

    return R;
}

void parseLine(std::string &line, double *R) {
    int j = 0;
    std::string s;
    for (auto d: line) {
        if (d == ' ') {
            R[j] = std::stod(s);
            s = "";
            ++j;
        } else s += d;
    }
    R[j] = std::stod(s);
}

double **load(std::fstream &file, int n) {
    double **R = creteSquareMatrix(n);
    for (int i = 0; i < n; ++i) {
        std::string line;
        std::getline(file, line);
        parseLine(line, R[i]);
    }
    return R;
}

int loadSize(std::fstream &file) {
    std::string s;
    std::getline(file, s);
    return std::stoi(s);
}

TEST_CASE("test 1") {
    std::fstream input("../lu_tests/test1/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test1/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 2") {
    std::fstream input("../lu_tests/test2/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test2/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 3") {
    std::fstream input("../lu_tests/test3/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test3/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 4") {
    std::fstream input("../lu_tests/test4/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test4/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 5") {
    std::fstream input("../lu_tests/test5/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test5/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 6") {
    std::fstream input("../lu_tests/test6/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test6/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 7") {
    std::fstream input("../lu_tests/test7/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test7/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 8") {
    std::fstream input("../lu_tests/test8/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test8/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 9") {
    std::fstream input("../lu_tests/test9/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test9/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}

TEST_CASE("test 10") {
    std::fstream input("../lu_tests/test10/input.txt", std::ios_base::in);
    std::fstream result("../lu_tests/test10/result.txt", std::ios_base::in);

    int n = loadSize(input);
    double **A = load(input, n);

    double **R = calculateInverseMatrix(A, n);

    double **expected = load(result, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            CHECK(expected[i][j] == round(R[i][j] * 10000) / 10000);
        }
    }

    input.close();
    result.close();

    delete[] A;
    delete[] R;
}