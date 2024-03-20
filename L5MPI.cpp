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
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);

	MPI_Status status;

	if (rank == 0)
	{
			
		double **matrix = new double*[msize];
		double *rightnum = new double[msize];

		double* truesolutions = GenerateSLAU(matrix, rightnum, msize);
		chrono::steady_clock::time_point start = chrono::steady_clock::now();

		int nocode = -1;
		int psize = 0;
		int remains = 0;
		int startindex = 0;
		double* mainstring;
		double* massiv;
		double* rightpart;

		for (int a = 0; a < msize - 1; ++a)
		{
			psize = (msize - 1 - a) / commsize;
			remains = (msize - 1 - a) % commsize;
			mainstring = PartOfMatrixToMassiv(matrix, msize, a, 1);

			for (int b = 1; b < commsize; ++b)
			{
				if (psize)
				{
					startindex = a + 1 + remains + (b * psize);
					massiv = PartOfMatrixToMassiv(matrix, msize, startindex, psize);

					MPI_Send(&a, 1, MPI_INT, b, 0, MPI_COMM_WORLD);
					MPI_Send(mainstring, msize, MPI_DOUBLE, b, 1, MPI_COMM_WORLD);
					MPI_Send(massiv, msize * psize, MPI_DOUBLE, b, 2, MPI_COMM_WORLD);
					//delete[] massiv;
				}
				else
				{
					MPI_Send(&nocode, 1, MPI_INT, b, 0, MPI_COMM_WORLD);
				}

			}

			//delete[] mainstring;
			double multiplier = 0;
		    for (int b = a + 1; b < (a + psize + remains + 1); ++b)
		    {
		        multiplier = matrix[b][a] / matrix[a][a];
		        for (int c = a; c < msize; ++c)
		        {
		       		matrix[b][c] += matrix[a][c] * multiplier * -1;
		        }
		        rightnum[b] += rightnum[a] * multiplier * -1;
		    }

			for (int b = (a + psize + remains + 1); b < msize; ++b)
			{
				multiplier = matrix[b][a] / matrix[a][a];
		        rightnum[b] += rightnum[a] * multiplier * -1;
			}

			if (psize)
			{
				for (int b = 1; b < commsize; ++b)
				{
					startindex = a + 1 + remains + (b * psize);
					MPI_Recv(massiv, msize * psize, MPI_DOUBLE, b, 0, MPI_COMM_WORLD, &status);

					MassivToPartOfMatrix(matrix, massiv, msize, startindex, psize);
					//delete[] massiv;
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

		chrono::steady_clock::time_point finish = chrono::steady_clock::now();
		chrono::milliseconds calculatetime = chrono::duration_cast<chrono::milliseconds>(finish - start);

		cout << "Time of running: " << calculatetime.count() << endl;
		if (TestSolves(truesolutions, solutions, msize)) cout << "True. " << endl;
		else cout << "False. ";
		
		for (int a = 0; a < msize; ++a) delete[] matrix[a];
		delete[] matrix;
		delete[] rightnum;
		delete[] truesolutions;
		delete[] solutions;
	}
	else
	{
		int index = 0;
		int kstrings = 0;
		int masslength = 0;
		double* mainstring;
		double* massiv;
		double** matrix;

		for (int a = 0; a < msize - 1; ++a)
		{
			MPI_Recv(&index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			if (index >= 0)
			{
				MPI_Probe(0, 2, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_DOUBLE, &masslength);

				mainstring = new double[msize];
				massiv = new double[masslength];
				MPI_Recv(mainstring, msize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(massiv, masslength, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);

				kstrings = masslength / msize;
				matrix = new double*[kstrings];
				for (int b = 0; b < kstrings; ++b)
				{
					matrix[b] = new double[msize];
				}

				MassivToPartOfMatrix(matrix, massiv, msize, 0, kstrings);
				
				double multiplier = 0;
				for (int b = 0; b < kstrings; ++b)
				{
					multiplier = matrix[b][index] / mainstring[index];
					for (int c = index; c < msize; ++c)
					{
						matrix[b][c] += mainstring[c] * multiplier * -1;
					}
				}
				
				delete[] mainstring;
				delete[] massiv;
				massiv = PartOfMatrixToMassiv(matrix, msize, 0, kstrings);
				for (int b = 0; b < kstrings; ++b) delete[] matrix[b];
				delete[] matrix;
				
				MPI_Send(massiv, masslength, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				delete[] massiv;
			}
		}	
	}
	
	MPI_Finalize();
	return 0;
}
