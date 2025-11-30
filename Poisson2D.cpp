#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <ctime>
#include <cmath>
#include <chrono>
#include <cassert>
#include <omp.h>
#include <mpi.h>

#include "CSRMatrix.h"
#include "Variant7.h"
#include "ConjugateGradient.h"
#include "MPINode.h"


size_t GetMilisecondsCount()
{
	using namespace std::chrono;
	return  duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// Rectangle split to N parts
std::vector<Domain> SplitDomain2D(int P, int& X_segments, int& Y_segments)
{
	if (P <= 0)
		throw std::invalid_argument("P should be positive even number");

	if (P == 1)
	{
		X_segments = Y_segments = 1;
		return std::vector<Domain>{{X_min, X_max, Y_min, Y_max}};
	}

	if (P % 2 == 1)
		throw std::invalid_argument("P should be positive even number");

	std::vector<Domain> subdomains;
	subdomains.reserve(P);

	int best_nx = 1, best_ny = P;
	double best_ratio = 10.0;

	// optimal n_x, n_y
	for (int i = 0; i <= std::log2(P); ++i)
	{
		int nx = 1 << i;
		int ny = P / nx;
		double ratio = static_cast<double>(nx) / ny;

		if (ratio >= 0.5 && ratio <= 2.0)
		{
			double deviation = std::fabs(std::log2(ratio));
			if (deviation < best_ratio)
			{
				best_ratio = deviation;
				best_nx = nx;
				best_ny = ny;
			}
		}
	}

	X_segments = best_nx;
	Y_segments = best_ny;

	double dx = (X_max - X_min) / best_nx;
	double dy = (Y_max - Y_min) / best_ny;

	for (int j = 0; j < best_ny; ++j)
	{
		for (int i = 0; i < best_nx; ++i)
		{
			Domain s;
			s.x_min = X_min + i * dx;
			s.x_max = X_min + (i + 1) * dx;
			s.y_min = Y_min + j * dy;
			s.y_max = Y_min + (j + 1) * dy;
			subdomains.push_back(s);
		}
	}

	return subdomains;
}

void OMPTest()
{
	int N, M; // X axis partitioned to M segments, Y - to N
	int NumThreads;

	std::cin >> M >> N >> NumThreads;
	std::cout << "OpenMP test with " << M << "x" << N << " grid and " << NumThreads << " threads" << std::endl;

	omp_set_num_threads(NumThreads);

	CSRMatrix A;
	std::vector<double> F, omega; // these are matrixes, just flatten
	CreateMatrixesV7(A, F, M + 1, N + 1);

	size_t avgTime = 0;
	int Passes = 3;
	for (int i = 0; i < Passes; i++)
	{
		size_t start = GetMilisecondsCount();
		omega = ConjugateGradient(A, F);
		size_t end = GetMilisecondsCount();

		std::cout << i << " run - " << end - start << " ms" << std::endl;

		avgTime += end - start;
	}
	avgTime /= Passes;

	std::cout << "average " << avgTime << " ms" << std::endl;
}

void MPITest(int argc, char** argv)
{
	int N, M, NumThreads; // X axis partitioned to M segments, Y - to N
	int X_segments, Y_segments; // Number of segments per axis
	int world_rank, world_size;
	MPINode node;
	// this is matrix, just flatten
	std::vector<double> local_omega;

	M = 6;
	N = 8;
	NumThreads = 1;

	//std::cin >> M >> N >> NumThreads;
	omp_set_num_threads(NumThreads);

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Rank of the process
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Total number of processes

	if (world_rank == 0)
	{
		std::cout << "MPI test" << std::endl;
		std::cout << world_size << " nodes with " << NumThreads << " threads in each" << std::endl;
		std::cout << M << "x" << N << " grid" << std::endl;
	}

	Domain InitialDomain{ X_min, X_max, Y_min, Y_max, M, N, M, N};
	node.CreateDomainInfo(InitialDomain, M, N);
	CreateMatrixesV7(node.A, node.F, node.m_Subdomain.Nx_total, node.m_Subdomain.Ny_total, node.m_Subdomain);
	
	size_t avgTime = 0;
	int Passes = 1;
	for (int i = 0; i < Passes; i++)
	{
		size_t start = GetMilisecondsCount();
		local_omega = node.ConjugateGradient();
		size_t end = GetMilisecondsCount();

		if (world_rank == 0)
			std::cout << i << " run - " << end - start << " ms" << std::endl;

		avgTime += end - start;
	}
	avgTime /= Passes;

	if (world_rank == 0)
		std::cout << "average " << avgTime << " ms" << std::endl;

	//MPINode::GatherOmega(local_omega, M, N, X_segments, Y_segments, true);

	//std::string ResultFileName = std::to_string(world_rank) + "Result" + std::to_string(Nx_local) + "x" + std::to_string(Ny_local) + ".txt";
	//PrintFlatMatrix(ResultFileName, local_omega, Ny_local, Nx_local);

	MPI_Finalize();
}

//extern int NumThreads;

int main(int argc, char** argv)
{
	//OMPTest();
	MPITest(argc, argv);
}
