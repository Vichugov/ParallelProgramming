#include<iostream>
#include<cstdlib>
#include<chrono>
#include<cmath>
#include<mpi.h>

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


double* PartOfMatrixToMassiv(double** matrix, int msize, int startindex, int kstrings)
{
	double* massiv = new double[msize * kstrings];
	int endindex = startindex + kstrings;
	for (int a = startindex; a < endindex; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			massiv[(a - startindex) * msize + b] = matrix[a][b];
		}
	}
	return massiv;
}


double* TakePartOfMassiv(double* massiv, int startindex, int kelements)
{
	double* part = new double[kelements];
	int endindex = startindex + kelements;
	for (int a = startindex; a < endindex; ++a)
	{
		part[a - startindex] = massiv[a];
	}
	return part;
}


void MassivToPartOfMatrix(double** matrix, double* massiv, int msize, int startindex, int kstrings)
{
	int endindex = startindex + kstrings;
	for (int a = startindex; a < endindex; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			matrix[a][b] = massiv[(a - startindex) * msize + b];
		}
	}
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
	int rank;
	int commsize;
	int msize = 1000;
	int fsize;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Status status;

	int psize = msize / commsize;
	int remain = msize % commsize;
	int startindex;
	int finishindex;
	double** matrix;
	double* rightnum;
	double* truesolutions;
	double** fragment;
	double* massiv;
	double* rightpart;
	chrono::steady_clock::time_point start;

	if (rank == 0)
	{
	
		matrix = new double*[msize];
		rightnum = new double[msize];

		truesolutions = GenerateSLAU(matrix, rightnum, msize);
		start = chrono::steady_clock::now();

		for (int a = 1; a < commsize; ++a)
		{
			massiv = PartOfMatrixToMassiv(matrix, msize, a * psize + remain, psize);
			rightpart = TakePartOfMassiv(rightnum, a * psize + remain, psize);
			MPI_Send(massiv, msize * psize, MPI_DOUBLE, a, 0, MPI_COMM_WORLD);
			MPI_Send(rightpart, psize, MPI_DOUBLE, a, 1, MPI_COMM_WORLD);
		}

		fsize = psize + remain;
		fragment = new double*[fsize];
		rightpart = new double[fsize];
		for (int a = 0; a < fsize; ++a)
		{
			fragment[a] = new double[msize];
			rightpart[a] = rightnum[a];
			for (int b = 0; b < msize; ++b)
			{
				fragment[a][b] = matrix[a][b];
			}
		}
		startindex = 0;
		finishindex = fsize;
	}
	else
	{
		fsize = psize;
		massiv = new double[msize * fsize];
		rightpart = new double[fsize];
		MPI_Recv(massiv, msize * fsize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(rightpart, fsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

		fragment = new double*[fsize];
		for (int a = 0; a < fsize; ++a)
		{
			fragment[a] = new double[msize];
		}
		MassivToPartOfMatrix(fragment, massiv, msize, 0, fsize);
		startindex = (psize + remain) + (psize * (rank - 1));
		finishindex = startindex + psize;
	}

	double* mainstring = new double[msize];
	double rightelement;
	int sender = 0;
	double multiplier;

	for (int a = 0; a < msize; ++a)
	{
		if (a >= startindex && a < finishindex)
		{
			for (int b = rank + 1; b < commsize; ++b)
			{
				MPI_Send(fragment[a - startindex], msize, MPI_DOUBLE, b, 2, MPI_COMM_WORLD);
				MPI_Send(&rightpart[a - startindex], 1, MPI_DOUBLE, b, 3, MPI_COMM_WORLD);
			}

			for (int b = a - startindex + 1; b < fsize; ++b)
			{
				multiplier = fragment[b][a] / fragment[a - startindex][a];
				for (int c = a; c < msize; ++c)
				{
					fragment[b][c] += fragment[a - startindex][c] * multiplier * -1;
				}
				rightpart[b] += rightpart[a - startindex] * multiplier * -1;
			}
		}
		else
		{
			if (a < startindex)
			{
				sender = (a - remain) / psize;
				MPI_Recv(mainstring, msize, MPI_DOUBLE, sender, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&rightelement, 1, MPI_DOUBLE, sender, 3, MPI_COMM_WORLD, &status);
				
				for (int b = 0; b < fsize; ++b)
				{
					multiplier = fragment[b][a] / mainstring[a];
					for (int c = a; c < msize; ++c)
					{
						fragment[b][c] += mainstring[c] * multiplier * -1;
					}
					rightpart[b] += rightelement * multiplier * -1;
				}
			}
		}
	}

	if (rank == 0)
	{
		for (int a = 0; a < fsize; ++a)
		{
			rightnum[a] = rightpart[a];
			for (int b = 0; b < msize; ++b)
			{
				matrix[a][b] = fragment[a][b];
			}
			delete[] fragment[a];
		}

		rightpart = new double[psize];

		for (int a = 1; a < commsize; ++a)
		{
			MPI_Recv(massiv, msize * psize, MPI_DOUBLE, a, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(rightpart, psize, MPI_DOUBLE, a, 1, MPI_COMM_WORLD, &status);
			MassivToPartOfMatrix(matrix, massiv, msize, fsize + (a - 1) * psize, psize);
			for (int b = 0; b < psize; ++b)
			{
				rightnum[fsize + (a - 1) * psize + b] = rightpart[b];
			}
		}

		double sum = 0; 
		double *solutions = new double[msize];
		solutions[msize - 1] = rightnum[msize - 1] / matrix[msize - 1][msize - 1];
		for (int a = msize - 2; a >= 0; --a)
		{
		    sum = rightnum[a];
		    for (int b = msize - 1; b > a; --b)
		    {
		        sum -= matrix[a][b] * solutions[b];
		    }
		    solutions[a] = sum / matrix[a][a];
		}

		chrono::steady_clock::time_point finish = chrono::steady_clock::now();
		chrono::milliseconds calculatetime = chrono::duration_cast<chrono::milliseconds>(finish - start);

		cout << "Time of running: " << calculatetime.count() << endl;
		if (TestSolves(truesolutions, solutions, msize)) cout << "True. " << endl;
		else cout << "False. " << endl;

		for (int a = 0; a < msize; ++a)
		{
			delete[] matrix[a];
		}

		delete[] matrix;
		delete[] rightnum;
		delete[] solutions;
	}
	else
	{
		massiv = PartOfMatrixToMassiv(fragment, msize, 0, psize);
		MPI_Send(massiv, msize * psize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(rightpart, psize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

		for (int a = 0; a < psize; ++a)
		{
			delete[] fragment[a];
		}
	}

	delete[] fragment;
	delete[] rightpart;
	delete[] massiv;
			
	MPI_Finalize();
	return 0;
}
