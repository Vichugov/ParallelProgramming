#include<iostream>
#include<cstdlib>

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

void CreateMatrix(int** amatrix, int** bmatrix, int msize)
{
	for (int a = 0; a < msize; ++a)
	{
		amatrix[a] = new int[msize];
		bmatrix[a] = new int[msize];
		for (int b = 0; b < msize; ++b)
		{
			amatrix[a][b] = rand() % 10;
			bmatrix[a][b] = rand() % 10;
		}
	}
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
		int bsize = msize / 2;

		int**** brmatrix = new int***[blocks];
		int**** bfmatrix = new int***[blocks];
		int**** bsmatrix = new int***[blocks];

		for (int a = 0; a < blocks; ++a)
		{
			brmatrix[a] = new int**[blocks];
			bfmatrix[a] = new int**[blocks];
			bsmatrix[a] = new int**[blocks];

			for (int b = 0; b < blocks; ++b)
			{
				brmatrix[a][b] = new int*[bsize];
				bfmatrix[a][b] = new int*[bsize];
				bsmatrix[a][b] = new int*[bsize];

				for (int c = 0; c < bsize; ++c)
				{
					brmatrix[a][b][c] = new int[bsize];
					bfmatrix[a][b][c] = new int[bsize];
					bsmatrix[a][b][c] = new int[bsize];
					for (int d = 0; d < bsize; ++d)
					{
						brmatrix[a][b][c][d] = 0;
						bfmatrix[a][b][c][d] = 0;
						bsmatrix[a][b][c][d] = 0;
					}
				}
				TakePartOfMatrix(rmatrix, brmatrix[a][b], msize, a * 2 + b);
				TakePartOfMatrix(fmatrix, bfmatrix[a][b], msize, a * 2 + b);
				TakePartOfMatrix(smatrix, bsmatrix[a][b], msize, a * 2 + b);
			}
		}
		
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
		
		for (int a = 0; a < blocks; ++a)
		{
			for (int b = 0; b < blocks; ++b)
			{
				for (int c = 0; c < bsize; ++c)
				{
					delete[] bfmatrix[a][b][c];
					delete[] bsmatrix[a][b][c];
					delete[] brmatrix[a][b][c];
				}
				delete[] bfmatrix[a][b];
				delete[] bsmatrix[a][b];
				delete[] brmatrix[a][b];
			}
			delete[] bfmatrix[a];
			delete[] bsmatrix[a];
			delete[] brmatrix[a];
		}
		delete[] bfmatrix;
		delete[] bsmatrix;
		delete[] brmatrix;
	}
	else
	{
		MultipleMatrixSizeTwo(fmatrix, smatrix, rmatrix);
	}
}


int main()
{
	int size = 4;
	int** amatrix = new int*[size];
	int** bmatrix = new int*[size];
	CreateMatrix(amatrix, bmatrix, size);
	ShowMatrix(amatrix, size);
	ShowMatrix(bmatrix, size);

	int** rmatrix = CreateEmptyMatrix(size);

	BlockStep(amatrix, bmatrix, rmatrix, size);
	//MultipleMatrix(amatrix, bmatrix, rmatrix, size);
	ShowMatrix(rmatrix, size);


	return 0;
}
