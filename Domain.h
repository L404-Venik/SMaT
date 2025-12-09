#pragma once
#include <vector>

struct Domain
{
	double x_min, x_max;
	double y_min, y_max;
	int Nx_local, Ny_local; // Nodes to calculate values in
	int Nx_total, Ny_total;	// Nodes including overlap with neighbour domains
};

namespace domain
{
	/**
	 * @brief Checks the possibility of a (M x N nodes) domain structured partition
	 * into P subdomains R x C satisfying the sides ratio constraint.
	 *
	 * @param P Total number of subdomains.
	 * @param M Segments number by X (M+1 nodes).
	 * @param N Segments number by Y (N+1 nodes).
	 * @param R the number of subdomains returned along the axis X.
	 * @param C the number of domains returned along the axis Y.
	 * @return true, if partition R x C found, false otherwise.
	 */
	bool FindOptimalPartitionRC(int P, int M, int N, int& R, int& C);

	/**
	 * @brief Splits the source domain into P subdomains with 2 nodes overlap.
	 *
	 * @param P Total number of subdomains.
	 * @param InitialDomain Information about domain to make partition of
	 * @return std::vector<Domain> of subdomains.
	 */
	std::vector<Domain> SplitDomain2D(int P, const Domain& InitialDomain);
}