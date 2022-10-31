
#include <iostream>
#include <fstream>



void inverseMatrix(double** matrix, int n)
{
	double** inverse = new double* [n];
	for (int i = 0; i < n; i++)
	{
		inverse[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			inverse[i][j] = 0;
		}
		inverse[i][i] = 1;
	}

	for (int i = 0; i < n; i++)
	{
		double temp = matrix[i][i];
		for (int j = 0; j < n; j++)
		{
			matrix[i][j] /= temp;
			inverse[i][j] /= temp;
		}

		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				temp = matrix[j][i];
				for (int k = 0; k < n; k++)
				{
					matrix[j][k] -= matrix[i][k] * temp;
					inverse[j][k] -= inverse[i][k] * temp;
				}
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrix[i][j] = inverse[i][j];
		}
	}

	for (int i = 0; i < n; i++)
	{
		delete[] inverse[i];
	}
	delete[] inverse;
}
int main() {
	std::ifstream file("input.txt");
	int n;
	file >> n;
	double** matrix = new double* [n];
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			file >> matrix[i][j];
		}
	}

	inverseMatrix(matrix, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::ofstream file2("output.txt");

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			file2 << matrix[i][j] << " ";
		}
		file2 << std::endl;
	}
	file2.close();
	file.close();
	for (int i = 0; i < n; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;

	return 0;
}










