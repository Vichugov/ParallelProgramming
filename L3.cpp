#include<iostream>
#include<cstdlib>

using namespace std;

int size;
int **matrix;

void ShowMatrix(double **matrix, double *rightnum)
{
    for (int a = 0; a < size; ++a)
    {
        for (int b = 0; b < size; ++b)
        {
            cout << matrix[a][b] << " ";
        }
        cout << " " << rightnum[a] << endl;
    }
    cout << endl;
}


void GenerateSLAU(double **matrix, double *rightnum)
{
    double *solutions = new double[size];
    
    for (int a = 0; a < size; ++a)
    {
        matrix[a] = new double[size];
        solutions[a] = rand() % 9 + 1;
    }

    for (int a = 0; a < size; ++a)
    {
        rightnum[a] = 0;

        for (int b = 0; b < size; ++b)
        {
            matrix[a][b] = rand() % 9 + 1;
            rightnum[a] += matrix[a][b] * solutions[b];
        }
    }

    ShowMatrix(matrix, solutions);
}


void SolveSLAU(double **matrix, double *rightnum)
{
    for (int a = 0; a < size - 1; ++a)
    {
        for (int b = a + 1; b < size; ++b)
        {
            double multiplier = (double)matrix[b][a] / (double)matrix[a][a];
            for (int c = a; c < size; ++c)
            {
                matrix[b][c] += matrix[a][c] * multiplier * -1;
            }
            rightnum[b] += rightnum[a] * multiplier * -1;
            ShowMatrix(matrix, rightnum);
        }
    }

    double *solutions = new double[size];
    solutions[size - 1] = rightnum[size - 1] / matrix[size - 1][size - 1];
    for (int a = size - 2; a >= 0; --a)
    {
        double sum = rightnum[a];
        for (int b = size - 1; b > a; --b)
        {
            sum -= matrix[a][b] * solutions[b];
        }
        solutions[a] = sum / matrix[a][a];
    }

    ShowMatrix(matrix, solutions);
}


int main()
{
    size = 10;
    double **matrix = new double*[size];
    double *rightnum = new double[size];

    GenerateSLAU(matrix, rightnum);
    ShowMatrix(matrix, rightnum);
    SolveSLAU(matrix, rightnum);
    ShowMatrix(matrix, rightnum);
}
