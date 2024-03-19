#include<iostream>
#include<cstdlib>
#include<chrono>
#include<cmath>

using namespace std;


void ShowMatrix(double **matrix, double *rightnum, int msize)
{
    for (int a = 0; a < msize; ++a)
    {
        for (int b = 0; b < msize; ++b)
        {
            cout << matrix[a][b] << " ";
        }
        cout << " " << rightnum[a] << endl;
    }
    cout << endl;
}


double* GenerateSLAU(double **matrix, double *rightnum, int msize)
{
    double *solutions = new double[msize];
    
    for (int a = 0; a < msize; ++a)
    {
        matrix[a] = new double[msize];
        solutions[a] = rand() % 9 + 1;
    }

    for (int a = 0; a < msize; ++a)
    {
        rightnum[a] = 0;

        for (int b = 0; b < msize; ++b)
        {
            matrix[a][b] = rand() % 9 + 1;
            rightnum[a] += matrix[a][b] * solutions[b];
        }
    }
	return solutions;
}


double* SolveSLAU(double **matrix, double *rightnum, int msize)
{

    for (int a = 0; a < msize - 1; ++a)
    {
    	{
    	    for (int b = a + 1; b < msize; ++b)
    	    {
    	        double multiplier = (double)matrix[b][a] / (double)matrix[a][a];
    	        int sum;
    	        for (int c = a; c < msize; ++c)
    	        {
    	       		matrix[b][c] += matrix[a][c] * multiplier * -1;
    	        }
    	        rightnum[b] += rightnum[a] * multiplier * -1;
    	    }
		}
    }

    double *solutions = new double[msize];
    solutions[msize - 1] = rightnum[msize - 1] / matrix[msize - 1][msize - 1];
    for (int a = msize - 2; a >= 0; --a)
    {
        double sum = rightnum[a];
        for (int b = msize - 1; b > a; --b)
        {
            sum -= matrix[a][b] * solutions[b];
        }
        solutions[a] = sum / matrix[a][a];
    }
	return solutions;
}


bool TestSolves(double* fsol, double* ssol, int msize)
{
	bool norma = true;
	for (int a = 0; a < msize; ++a)
	{
		if (abs((double)fsol[a] - (double)ssol[a]) > 0.001)
		{
			norma = false;
		}
	}
	return norma;
}

int main()
{
    int size = 1000;
    double **matrix = new double*[size];
    double *rightnum = new double[size];

    double* truesolutions = GenerateSLAU(matrix, rightnum, size);

	chrono::steady_clock::time_point start = chrono::steady_clock::now();
    double* solutions = SolveSLAU(matrix, rightnum, size);
	chrono::steady_clock::time_point finish = chrono::steady_clock::now();
	chrono::milliseconds calculatetime = chrono::duration_cast<chrono::milliseconds>(finish - start);

	cout << "Time of running: " << calculatetime.count() << endl;
	if (TestSolves(truesolutions, solutions, size)) cout << "True. " << endl;
	else cout << "False. ";
	
	for (int a = 0; a < size; ++a) delete[] matrix[a];
	delete[] matrix;
	delete[] rightnum;
	delete[] truesolutions;
	delete[] solutions;
}
