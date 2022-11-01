#ifndef MATRIX_HPP
#define MATIX_HPP
#include <iostream>
using namespace std;
template <typename T>
class Matrix
{
private:
    int size;
    T **matr;

public:
    Matrix(T **matr, int i, int j);
    Matrix(const std::string &path);
    Matrix(const Matrix &m);
    bool INV();
    double DET(int n);
    void getCfactor(int p, int q, int n);

    ~Matrix();
};
#endif

// TODO
template <typename T>
Matrix<T>::Matrix(T **matr, int i, int j)
{
}

// TODO
template <typename T>
Matrix<T>::Matrix(const Matrix &m)
{
}

// TODO
template <typename T>
Matrix<T>::Matrix(const std::string &path)
{
}

#define N 5

template <typename T>
void Matrix<T>::getCfactor(int p, int q, int n)
{
    T **t(this);
    int i = 0, j = 0;
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++) // Copy only those elements which are not in given row r and column c:
        {
            if (r != p && c != q)
            {
                t[i][j++] = M[r][c]; // If row is filled increase r index and reset c index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
template <typename T>
double Matrix<T>::DET(int n) // to find determinant
{
    double D = 0;
    if (n == 1)
        return this->matr;
    T t[this->size][this->size]; // store cofactors
    int s = 1;                   // store sign multiplier
    // To Iterate each element of first row
    for (int f = 0; f < n; f++)
    {
        // For Getting Cofactor of M[0][f] do
        getCfactor(0, f, n);
        D += s * M[0][f] * DET(n - 1);
        s = -s;
    }
    return D;
}
void ADJ(int M[N][N], int adj[N][N])
// to find adjoint matrix
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }
    int s = 1,
        t[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // To get cofactor of M[i][j]
            getCfactor(M, t, i, j, N);
            s = ((i + j) % 2 == 0) ? 1 : -1;   // sign of adj[j][i] positive if sum of row and column indexes is even.
            adj[j][i] = (s) * (DET(t, N - 1)); // Interchange rows and columns to get the transpose of the cofactor matrix
        }
    }
}
template <typename T>
bool Matrix<T>::INV()
{
    int det = DET(M, N);
    if (det == 0)
    {
        cout << "can't find its inverse";
        return false;
    }
    int adj[N][N];
    ADJ(M, adj);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inv[i][j] = adj[i][j] / float(det);
    return true;
}
template <class T>
void print(T A[N][N]) // print the matrix.
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}