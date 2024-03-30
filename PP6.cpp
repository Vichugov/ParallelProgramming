#include<iostream>
#include<cstdlib>
#include<chrono>
#include<cmath>
#include<tbb/tbb.h>

using namespace tbb;


class PForSLAU
{
	double** matrix;
	double* rightnum;
	int msize;
	int a;
public:

	PForSLAU(double** pmatrix, double* prightnum, int pmsize, int pa) : matrix(pmatrix), rightnum(prightnum), msize(pmsize), a(pa)
	{

	}

	void operator()(const blocked_range<int>& range) const
	{
		for (int b = range.begin(); b < range.end(); ++b)
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
};


void ShowMatrix(double** matrix, double* rightnum, int msize)
{
	for (int a = 0; a < msize; ++a)
	{
		for (int b = 0; b < msize; ++b)
		{
			std::cout << matrix[a][b] << " ";
		}
		std::cout << " " << rightnum[a] << std::endl;
	}
	std::cout << std::endl;
}


double* GenerateSLAU(double** matrix, double* rightnum, int msize)
{
	double* solutions = new double[msize];

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


double* SolveSLAU(double** matrix, double* rightnum, int msize)
{
	for (int a = 0; a < msize - 1; ++a)
	{
		parallel_for(blocked_range<int>(a + 1, msize), PForSLAU(matrix, rightnum, msize, a));
	}

	double* solutions = new double[msize];
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
	int size = 2000;
	double** matrix = new double* [size];
	double* rightnum = new double[size];

	double* truesolutions = GenerateSLAU(matrix, rightnum, size);

	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	double* solutions = SolveSLAU(matrix, rightnum, size);
	std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
	std::chrono::milliseconds calculatetime = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

	std::cout << "Time of running: " << calculatetime.count() << std::endl;
	if (TestSolves(truesolutions, solutions, size)) std::cout << "True. " << std::endl;
	else std::cout << "False. ";

	for (int a = 0; a < size; ++a) delete[] matrix[a];
	delete[] matrix;
	delete[] rightnum;
	delete[] truesolutions;
	delete[] solutions;

	char z;
	std::cin >> z;

	return 0;
}