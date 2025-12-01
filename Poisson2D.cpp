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

	//std::string ResultFileName = "Result" + std::to_string(M) + "x" + std::to_string(N) + ".txt";
	//PrintFlatMatrix(ResultFileName, omega, N + 1, M + 1);
}

void MPITest(int argc, char** argv)
{
	int N, M, NumThreads; // X axis partitioned to M segments, Y - to N
	int world_rank, world_size;
	MPINode node;
	// this is matrix, just flatten
	std::vector<double> local_omega;

	M = 400;
	N = 600;
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

	Domain InitialDomain{ X_min, X_max, Y_min, Y_max, M + 1, N + 1, M + 1, N + 1};
	node.CreateDomainInfo(InitialDomain);
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

	//std::string ResultFileName = std::to_string(world_rank) + "Result.txt";
	//PrintFlatMatrix(ResultFileName, local_omega, node.m_Subdomain.Ny_total, node.m_Subdomain.Nx_total);

	MPI_Finalize();
}

int main(int argc, char** argv)
{
	//OMPTest();
	MPITest(argc, argv);
}
