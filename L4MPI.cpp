#include<iostream>
#include<cstdlib>
#include<chrono>
#include<cmath>
#include<mpi.h>

using namespace std;


void ShowMatrix(int** matrix, int msize)
{
	for (int a = 0; a < msize; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			cout << matrix[a][b] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


int** CreateEmptyMatrix(int msize)
{
	int** matrix = new int*[msize];
	for (int a = 0; a < msize; ++a)
	{
		matrix[a] = new int[msize];
		for (int b = 0; b < msize; ++b)
		{
			matrix[a][b] = 0;
		}
	}
	return matrix;
}

int** CreateRandomMatrix(int msize)
{
	int** matrix = new int*[msize];
	for (int a = 0; a < msize; ++a)
	{
		matrix[a] = new int[msize];
		for (int b = 0; b < msize; ++b)
		{
			matrix[a][b] = rand() % 10;
		}
	}
	return matrix;
}


void TakePartOfMatrix(int** matrix, int** pmatrix, int msize, int part)
{
	int psize = msize / 2;
	for (int a = 0; a < psize; ++a) pmatrix[a] = new int[psize];
	switch(part)
	{
		case 0:
			for (int a = 0; a < psize; ++a)
			{
				for (int b = 0; b < psize; ++b) pmatrix[a][b] = matrix[a][b];
			}
		break;
		case 1:
			for (int a = 0; a < psize; ++a)
			{
				for (int b = psize; b < msize; ++b) pmatrix[a][b - psize] = matrix[a][b];
			}
		break;
		case 2:
			for (int a = psize; a < msize; ++a)
			{
				for (int b = 0; b < psize; ++b) pmatrix[a - psize][b] = matrix[a][b];
			}
		break;
		case 3:
			for (int a = psize; a < msize; ++a)
			{
				for (int b = psize; b < msize; ++b) pmatrix[a - psize][b - psize] = matrix[a][b];
			}
		break;
	}
}


void InsertPartIntoMatrix(int** matrix, int** pmatrix, int msize, int part)
{
	int psize = msize / 2;
	switch(part)
	{
		case 0:
			for (int a = 0; a < psize; ++a)
			{
				for (int b = 0; b < psize; ++b) matrix[a][b] = pmatrix[a][b];
			}
		break;
		case 1:
			for (int a = 0; a < psize; ++a)
			{
				for (int b = psize; b < msize; ++b) matrix[a][b] = pmatrix[a][b - psize];
			}
		break;
		case 2:
			for (int a = psize; a < msize; ++a)
			{
				for (int b = 0; b < psize; ++b) matrix[a][b] = pmatrix[a - psize][b];
			}
		break;
		case 3:
			for (int a = psize; a < msize; ++a)
			{
				for (int b = psize; b < msize; ++b) matrix[a][b] = pmatrix[a - psize][b - psize];
			}
		break;
	}
}


int**** GetBlocks(int** matrix, int msize, int blocks)
{
	int bsize = msize / blocks;
	int**** bmatrix = new int***[blocks];

	for (int a = 0; a < blocks; ++a)
	{
		bmatrix[a] = new int**[blocks];
		for (int b = 0; b < blocks; ++b)
		{
			bmatrix[a][b] = new int*[bsize];
			TakePartOfMatrix(matrix, bmatrix[a][b], msize, a * 2 + b);
		}
	}
	
	return bmatrix;
}


void DelBlocks(int**** bmatrix, int msize, int blocks)
{
	int bsize = msize / 2;
	for (int a = 0; a < blocks; ++a)
	{
		for (int b = 0; b < blocks; ++b)
		{
			for (int c = 0; c < bsize; ++c)
			{
				delete[] bmatrix[a][b][c];
			}
			delete[] bmatrix[a][b];
		}
		delete[] bmatrix[a];
	}
	delete[] bmatrix;
}


void AddMatrix(int** rmatrix, int** smatrix, int msize)
{
	for (int a = 0; a < msize; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			rmatrix[a][b] += smatrix[a][b];
		}
	}
}


void MultipleMatrixSizeTwo(int** fmatrix, int** smatrix, int** rmatrix)
{
	rmatrix[0][0] = fmatrix[0][0] * smatrix[0][0] + fmatrix[0][1] * smatrix[1][0];
	rmatrix[0][1] = fmatrix[0][0] * smatrix[0][1] + fmatrix[0][1] * smatrix[1][1];
	rmatrix[1][0] = fmatrix[1][0] * smatrix[0][0] + fmatrix[1][1] * smatrix[1][0];
	rmatrix[1][1] = fmatrix[1][0] * smatrix[0][1] + fmatrix[1][1] * smatrix[1][1];
}


void MultipleMatrix(int** fmatrix, int** smatrix, int** rmatrix, int msize)
{
	for (int a = 0; a < msize; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			rmatrix[a][b] = 0;
			for (int c = 0; c < msize; ++c)
			{
				rmatrix[a][b] += fmatrix[a][c] * smatrix[c][b];
			}
		}
	}
}


void BlockStep(int** fmatrix, int** smatrix, int** rmatrix, int msize)
{
	if (msize > 2)
	{
		int blocks = 2;
		int bsize = msize / blocks;

		int**** bfmatrix = GetBlocks(fmatrix, msize, blocks);
		int**** bsmatrix = GetBlocks(smatrix, msize, blocks);
		int**** brmatrix = GetBlocks(rmatrix, msize, blocks);

		int** multipleresult = CreateEmptyMatrix(bsize);

		for (int a = 0; a < blocks; ++a)
		{
			for (int b = 0; b < blocks; ++b)
			{
				for (int c = 0; c < blocks; ++c)
				{
					BlockStep(bfmatrix[a][c], bsmatrix[c][b], multipleresult, bsize);
					AddMatrix(brmatrix[a][b], multipleresult, bsize);
					for (int d = 0; d < bsize; ++d) 
					{
						for (int e = 0; e < bsize; ++e) multipleresult[d][e] = 0;
					}
				}
				InsertPartIntoMatrix(rmatrix, brmatrix[a][b], msize, a * 2 + b);
			}
		}
		DelBlocks(bfmatrix, msize, blocks);
		DelBlocks(bsmatrix, msize, blocks);
		DelBlocks(brmatrix, msize, blocks);
		for (int a = 0; a < bsize; ++a) 
		{
			delete[] multipleresult[a];
		}
		delete[] multipleresult;
	}
	else
	{
		MultipleMatrixSizeTwo(fmatrix, smatrix, rmatrix);
	}
}


bool TestResult(int** fresult, int** sresult, int rsize)
{
	bool norma = true;
	for (int a = 0; a < rsize; ++a)
	{
		for (int b = 0; b < rsize; ++b)
		{
			if (fresult[a][b] != sresult[a][b])
			{
				norma = false;
			}
		}
	}
	return norma;
}


int* MatrixToMassiv(int** matrix, int msize)
{
	int* massiv = new int[msize * msize];
	for (int a = 0; a < msize; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			massiv[a * msize + b] = matrix[a][b];
		}
	}
	return massiv;
}


int** MassivToMatrix(int* massiv, int msize)
{
	int** matrix = new int*[msize];
	for (int a = 0; a < msize; ++a)
	{
		matrix[a] = new int[msize];
		for (int b = 0; b < msize; ++b)
		{
			matrix[a][b] = massiv[a * msize + b];
		}
	}
	return matrix;
}

int main()
{
	int blocks = 2;
	int rank;
	int commsize;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);

	int nomer;
	MPI_Status status;

	if (rank == 0)
	{
		int size = 256;
		int bsize = size / blocks;

		int** fmatrix = CreateRandomMatrix(size);
		int** smatrix = CreateRandomMatrix(size);

		int** tmatrix = CreateEmptyMatrix(size);
		int** rmatrix = CreateEmptyMatrix(size);

		chrono::steady_clock::time_point start = chrono::steady_clock::now();
		int**** bfmatrix = GetBlocks(fmatrix, size, blocks);
		int**** bsmatrix = GetBlocks(smatrix, size, blocks);

		for (int a = 0; a < blocks; ++a)
		{
			for (int b = 0; b < blocks; ++b)
			{
				if (a + b)
				{
					int* fmassiv = MatrixToMassiv(fmatrix, size);
					int* smassiv = MatrixToMassiv(smatrix, size);

					MPI_Send(fmassiv, size * size, MPI_INT, a * 2 + b, 0, MPI_COMM_WORLD);
					MPI_Send(smassiv, size * size, MPI_INT, a * 2 + b, 1, MPI_COMM_WORLD);
				}
			}
		}

		int** tempmatrix = CreateEmptyMatrix(bsize);
		int** nnmatrix = CreateEmptyMatrix(bsize);
		for (int a = 0; a < blocks; ++a)
		{
			BlockStep(bfmatrix[0][a], bsmatrix[a][0], tempmatrix, bsize);
			AddMatrix(nnmatrix, tempmatrix, bsize);
			for (int b = 0; b < bsize; ++b)
			{
				for (int c = 0; c < bsize; ++c)
				{
					tempmatrix[b][c] = 0;
				}
			}
		}
		InsertPartIntoMatrix(rmatrix, nnmatrix, size, 0);

		int* massiv = new int[bsize * bsize];
		for (int a = 0; a < blocks; ++a)
		{
			for (int b = 0; b < blocks; ++b)
			{
				if (a + b)
				{
					MPI_Recv(massiv, bsize * bsize, MPI_INT, a * 2 + b, 2, MPI_COMM_WORLD, &status);
					delete[] tempmatrix;
					tempmatrix = MassivToMatrix(massiv, bsize);
					InsertPartIntoMatrix(rmatrix, tempmatrix, size, a * 2 + b);
				}
			}
		}

		chrono::steady_clock::time_point finish = chrono::steady_clock::now();
		chrono::milliseconds calculatetime = chrono::duration_cast<chrono::milliseconds>(finish - start);
		MultipleMatrix(fmatrix, smatrix, tmatrix, size);
		
		cout << "Time of running: " << calculatetime.count() << endl;
		cout << TestResult(rmatrix, tmatrix, size) << endl;
	}
	else
	{
		int** fmatrix;
		int** smatrix;
		int* massiv;

		int masssize = 0;

		MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &masssize);
		int msize = sqrt(masssize);
		massiv = new int[masssize];

		MPI_Recv(massiv, masssize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		fmatrix = MassivToMatrix(massiv, msize);

		MPI_Recv(massiv, masssize, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		smatrix = MassivToMatrix(massiv, msize);

		int**** fbmatrix = GetBlocks(fmatrix, msize, blocks);
		int**** sbmatrix = GetBlocks(smatrix, msize, blocks);

		int bsize = msize / blocks;
		int** rblock = CreateEmptyMatrix(bsize);
		int** tempmatrix = CreateEmptyMatrix(bsize);

		int hor = rank / 2;
		int ver = rank % 2;

		for (int a = 0; a < blocks; ++a)
		{
			BlockStep(fbmatrix[hor][a], sbmatrix[a][ver], tempmatrix, bsize);
			AddMatrix(rblock, tempmatrix, bsize);
			for (int b = 0; b < bsize; ++b)
			{
				for (int c = 0; c < bsize; ++c)
				{
					tempmatrix[b][c] = 0;
				}
			}
		}

		delete[] massiv;
		massiv = MatrixToMassiv(rblock, bsize);
		MPI_Send(massiv, bsize * bsize, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
	return 0;
}
