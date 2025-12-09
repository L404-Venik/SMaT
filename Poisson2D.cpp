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

struct TestParameters
{
	int M, N, NumThreads, Passes;
};

void OMPTest(int argc, char** argv, const TestParameters& Param)
{
	std::cout << "OpenMP test with " << omp_get_max_threads() << " threads" << std::endl;
	std::cout << Param.M << "x" << Param.N << " grid" << std::endl;

	CSRMatrix A;
	std::vector<double> F, omega; // these are matrixes, just flatten
	CreateMatrixesV7(A, F, Param.M + 1, Param.N + 1);

	size_t avgTime = 0;
	for (int i = 0; i < Param.Passes; i++)
	{
		size_t start = GetMilisecondsCount();
		omega = ConjugateGradient(A, F);
		size_t end = GetMilisecondsCount();

		std::cout << i << " run - " << end - start << " ms" << std::endl;

		avgTime += end - start;
	}
	avgTime /= Param.Passes;

	std::cout << "average " << avgTime << " ms" << std::endl;

	//std::string ResultFileName = "Result.txt";
	//PrintFlatMatrix(ResultFileName, omega, Param.M + 1, Param.N + 1);
}

void MPITest(int argc, char** argv, const TestParameters& Param)
{
	int world_rank, world_size;
	MPINode node;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Rank of the process
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Total number of processes

	if (world_rank == 0)
	{
		std::cout << "MPI test" << std::endl;
		std::cout << world_size << " nodes with " << omp_get_max_threads() << " threads in each" << std::endl;
		std::cout << Param.M << "x" << Param.N << " grid" << std::endl;
	}

	Domain InitialDomain{ X_min, X_max, Y_min, Y_max, Param.M + 1, Param.N + 1, Param.M + 1, Param.N + 1 };
	if (!node.CreateDomainInfo(InitialDomain))
	{
		std::cout << "Couldn't create node domain info. Finishing programm" << std::endl;
		MPI_Finalize();
		return;
	}
	CreateMatrixesV7(node.A, node.F, node.m_Subdomain.Nx_total, node.m_Subdomain.Ny_total, node.m_Subdomain);

	size_t avgTime = 0;
	for (int i = 0; i < Param.Passes; i++)
	{
		size_t start = GetMilisecondsCount();
		int iter = node.ConjugateGradient();
		size_t end = GetMilisecondsCount();

		if (world_rank == 0)
		{
			if (i == 0 && iter > 0)
				std::cout << "Converged in " << iter << " iterations" << std::endl;

			std::cout << i << " run - " << end - start << " ms" << std::endl;
		}

		avgTime += end - start;
	}
	avgTime /= Param.Passes;

	if (world_rank == 0)
		std::cout << "average " << avgTime << " ms" << std::endl;

	//MPINode::GatherOmega(local_omega, Param.M, Param.N, X_segments, Y_segments, true);

	MPI_Finalize();
}

int main(int argc, char** argv)
{
	TestParameters Param;

	// for local testing
	std::cin >> Param.M >> Param.N >> Param.NumThreads;


	// for testing on Polus
	/*if (argc < 3) 
	{
		std::cerr << "Usage: program <int1> <int2> <int3>\n";
		return -1;
	}
	Param.M = std::stoi(argv[1]);
	Param.N = std::stoi(argv[2]);
	Param.NumThreads = std::atoi(argv[3]);*/

	Param.Passes = 5;
	omp_set_num_threads(Param.NumThreads);

	OMPTest(argc, argv, Param);
	//MPITest(argc, argv, Param);

	return 0;
}
