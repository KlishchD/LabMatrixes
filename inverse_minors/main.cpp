#include <iostream>
#include <fstream>
using namespace std;
double **init_matrix(const int &size)
{
    double **M = new double *[size];
    for (int i = 0; i < size; i++)
    {
        M[i] = new double[size];
    }
    return M;
}
void delete_matrix(double **M, const int &size)
{
    for (int i = 0; i < size; i++)
    {
        delete[] M[i];
    }
    delete[] M;
}
void getCfactor(double **M, const int &size, double **t, int p, int q, int n)
{
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
int DET(double **M, const int &size, int n) // to find determinant
{
    int D = 0;
    if (n == 1)
        return M[0][0];
    double **t = init_matrix(size); // store cofactors
    int s = 1;                      // store sign multiplier //
    // To Iterate each element of first row
    for (int f = 0; f < n; f++)
    {
        // For Getting Cofactor of M[0][f] do
        getCfactor(M, size, t, 0, f, n);
        D += s * M[0][f] * DET(t, size, n - 1);
        s = -s;
    }
    delete_matrix(t, size);
    return D;
}
void ADJ(double **M, double **adj, const int &size)
// to find adjoint matrix
{
    if (size == 1)
    {
        adj[0][0] = 1;
        return;
    }
    int s = 1;
    double **t = init_matrix(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // To get cofactor of M[i][j]
            getCfactor(M, size, t, i, j, size);
            s = ((i + j) % 2 == 0) ? 1 : -1;            // sign of adj[j][i] positive if sum of row and column indexes is even.
            adj[j][i] = (s) * (DET(t, size, size - 1)); // Interchange rows and columns to get the transpose of the cofactor matrix
        }
    }
    delete_matrix(t, size);
}
double **INV(double **M, const int &size)
{
    int det = DET(M, size, size);
    if (det == 0)
    {
        return nullptr;
    }
    double **inv = init_matrix(size);
    double **adj = init_matrix(size);
    ADJ(M, adj, size);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            inv[i][j] = adj[i][j] / float(det);
    return inv;
}
template <class T>
void print(T **A, int size) // print the matrix.
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}
double **read_from_file(int &size, const std::string &path)
{
    fstream file(path);
    if (!file)
    {
        return nullptr;
    }
    file >> size;
    double **M = init_matrix(size);
    double tmp;
    int i = 0, j = 0;
    while (file >> tmp)
    {
        // *(M + i * size + j) = tmp;
        M[i][j] = tmp;
        j++;
        if (j >= size)
        {
            j = 0;
            i++;
        }
    }
    file.close();
    return M;
}
void print_matrix(double **M, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << M[i][j] << ' ';
        }
        cout << endl;
    }
}
bool is_equal(double **M, double **T, const int &size, const string &path)
{
    if (!M && !T)
    {
        return true;
    }
    else if (!M || !T)
    {
        return false;
    }
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            // cout << M[i][j] << ', ' << T[i][j] << endl;
            if (M[i][j] != T[i][j])
            {
                string err = "Inverse matrix is not correct " + path;
                throw std::invalid_argument(err);
            }
        }

    return true;
}
int main()
{

    const string tests[] = {"test1", "test2", "test3", "test4", "test5"};
    for (auto path : tests)
    {
        // read the problem and find the invers
        int size;
        double **M = read_from_file(size, path + "/problem.txt"), **inv = INV(M, size);
        // if inverse exists check against answer.txt
        double **answer = read_from_file(size, path + "/answer.txt");
        try
        {
            is_equal(answer, inv, size, path);
            cout << "Inverse matrix is correct " + path << endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
        if (inv)
        {
            delete_matrix(answer, size);
            delete_matrix(inv, size);
        }
        delete_matrix(M, size);
    }
    return 0;
}