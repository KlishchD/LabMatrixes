//
// Created by Dmytro Klishch on 10/31/22.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <fstream>
#include <iostream>
#include "doctest.h"

void calculateLNonFirstColumn(int j, double **A, double **L, double **U, int n, bool &error) {
    for (int i = j; i < n; ++i) {
        L[i][j] = A[i][j];
        for (int k = 0; k < j; ++k) {
            L[i][j] -= U[k][j] * L[i][k];
        }
    }
}

void calculateLFirstColumn(double **A, double **L, int n, bool &error) {
    for (int i = 0; i < n; ++i) {
        L[i][0] = A[i][0];
    }
}

void calculateLColumn(int j, double **A, double **L, double **U, int n, bool &error) {
    if (j == 0) calculateLFirstColumn(A, L, n, error);
    else calculateLNonFirstColumn(j, A, L, U, n, error);
}

void calculateUNonFirstRow(int i, double **A, double **L, double **U, int n, bool &error) {
    for (int j = i; j < n; ++j) {
        if (L[i][i] == 0) {
            error = 1;
            return;
        }
        U[i][j] = A[i][j];
        for (int k = 0; k < i; ++k) {
            U[i][j] -= L[i][k] * U[k][j];
        }
        U[i][j] /= L[i][i];
    }
}

void calculateUFirstRow(double **A, double **U, int n, bool &error) {
    if (A[0][0] == 0) {
        error = 1;
        return;
    }
    for (int j = 0; j < n; ++j) {
        U[0][j] = A[0][j] / A[0][0];
    }
}

void calculateURow(int i, double **A, double **L, double **U, int n, bool &error) {
    if (i == 0) calculateUFirstRow(A, U, n, error);
    else calculateUNonFirstRow(i, A, L, U, n, error);
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

double **calculateInverseMatrix(double **A, int n, bool &error) {
    double **L = creteSquareMatrix(n);
    double **U = creteSquareMatrix(n);
    double **R = creteSquareMatrix(n);

    for (int i = 0; i < n && !error; ++i) {
        calculateLColumn(i, A, L, U, n, error);
        calculateURow(i, A, L, U, n, error);
    }

    for (int i = 0; i < n && !error; ++i) {
        calculateRColumn(i, R, L, U, n);
    }

    delete[] L;
    delete[] U;

    if (error) {
        delete[] R;
        return nullptr;
    }

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
    std::cout << "::" << s << std::endl;
    return std::stoi(s);
}

double **multiply(double **A, double **B, int n) {
    double **R = creteSquareMatrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                R[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return R;
}

double EPS = 1e-9;

TEST_CASE("legitimate matrices") {
    for (int test = 0; test <= 10; ++test) {
        std::string caseName = "file_" + std::to_string(test);
        SUBCASE(caseName.c_str()) {
            std::fstream input("lu/tests/test" + std::to_string(test) + ".txt", std::ios_base::in);

            int n = loadSize(input);
            double **A = load(input, n);

            bool error = 0;

            double **R = calculateInverseMatrix(A, n, error);

            double **multiplied = multiply(A, R, n);

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::cout << multiplied[i][j] << " ";
                }
                std::cout << std::endl;
            }

            for (int i = 0; i < n; ++i) {
                CHECK(std::abs(multiplied[i][i] - 1) < EPS);
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    CHECK(multiplied[i][j] < EPS);
                }
            }

            input.close();

            delete[] A;
            delete[] R;
        }
    }
}

TEST_CASE("non-legitimate matrices") {
    for (int test = 11; test <= 13; ++test) {
        std::string caseName = "file_" + std::to_string(test);
        SUBCASE(caseName.c_str()) {
            std::fstream input("lu/tests/test" + std::to_string(test) + ".txt", std::ios_base::in);

            int n = loadSize(input);
            double **A = load(input, n);

            bool error = 0;

            double **R = calculateInverseMatrix(A, n, error);

            CHECK(error);

            input.close();

            delete[] A;
            delete[] R;
        }
    }
}
