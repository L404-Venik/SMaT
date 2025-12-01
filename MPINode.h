#pragma once
#include <algorithm>
#include <array>
#include "CSRMatrix.h"
#include "Domain.h"

enum class Direction { LEFT, RIGHT, UP, DOWN };

struct NeighborInfo
{
	int NeighborRank;            // which MPI rank it communicates with
	Direction direction;
	std::vector<int> sendIdx;    // local indices you send
	std::vector<int> recvIdx;    // local ghost indices to fill
	
	std::vector<double> tmpSend, tmpRecv;
};

class MPINode
{
	void ExchangeGhosts(std::vector<double>& v);
	double DotProduct(const std::vector<double>& x, const std::vector<double>& y);

	void BuildNeighborInfo(int world_rank, int X_segments, int Y_segments);
	int GetHorizontalNeighboursCount();
	int GetVertialNeighboursCount();
	int GetNeighboursCount();

	std::vector<NeighborInfo> m_vNeighbors;
	std::array<bool,4> m_aHasNeighbour;

	int m_iMaxIter;
public:

	CSRMatrix A;
	std::vector<double> F;
	Domain m_Subdomain;

	void CreateDomainInfo(const Domain& InitialDomain);

	std::vector<double> ConjugateGradient();
	static void GatherOmega(const std::vector<double>& omega, int M, int N, int X_segments, int Y_segments, bool bSave = false);
};