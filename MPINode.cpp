#include <assert.h>
#include <string>

#include "MPINode.h"
#include "mpi.h"

NeighborDirection opposite(NeighborDirection d)
{
	switch (d)
	{
	case NeighborDirection::LEFT:  return NeighborDirection::RIGHT;
	case NeighborDirection::RIGHT: return NeighborDirection::LEFT;
	case NeighborDirection::UP:    return NeighborDirection::DOWN;
	case NeighborDirection::DOWN:  return NeighborDirection::UP;
	}
}

void MPINode::ExchangeGhosts(std::vector<double>& v)
{
	std::vector<MPI_Request> requests;
	requests.reserve(m_vNeighbors.size() * 2);

	for (auto& nb : m_vNeighbors)
	{
		int count = nb.sendIdx.size();

		// Pack send buffer
		nb.tmpRecv.resize(count);
		nb.tmpSend.resize(count);

		for (int i = 0; i < count; i++)
			nb.tmpSend[i] = v[nb.sendIdx[i]];

		MPI_Request req1, req2;

		MPI_Irecv(
			nb.tmpRecv.data(), count,
			MPI_DOUBLE, nb.NeighborRank, 100 + (int)nb.direction,
			MPI_COMM_WORLD, &req1);

		MPI_Isend(
			nb.tmpSend.data(), count,
			MPI_DOUBLE, nb.NeighborRank, 100 + (int)opposite(nb.direction),
			MPI_COMM_WORLD, &req2);

		// We need recvBuf after wait, so capture by lambda
		requests.push_back(req1);
		requests.push_back(req2);
	}

	// Wait for all to finish
	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

	// Copy received ghosts into local vector
	for (auto& nb : m_vNeighbors)
	{
		for (size_t i = 0; i < nb.recvIdx.size(); i++)
			v[nb.recvIdx[i]] = nb.tmpRecv[i];

		nb.tmpSend.clear();
		nb.tmpRecv.clear();
	}
}


std::vector<double> MPINode::ConjugateGradient()
{
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int n = A.m_iRows;
	int M = std::sqrt(n);
	const int max_iter = n * world_size;
	const double delta = 0.01;

	std::vector<double> omega(n, 0.0);
	std::vector<double> r = F;          // r0 = F - A*x = F
	std::vector<double> p;
	std::vector<double> Ap(n, 0.0), z(n, 0.0);
	std::vector<double> D = A.GetDiagonal();

	for (int i = 0; i < n; ++i)
		z[i] = r[i] / D[i];

	p = z;

	if (world_rank == 0)
	{
		std::ofstream os("residual.txt", std::ios::out);
		os.close();
	}


	double rz_local = DotProduct(z, r);
	double rz_old = 0.0;
	MPI_Allreduce(&rz_local, &rz_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int it = 0; it < max_iter; it++)
	{
		ExchangeGhosts(p);

		Ap = A.VectorMultiply(p);


		double pAp_local = DotProduct(p, Ap);
		double pAp = 0.0;
		MPI_Allreduce(&pAp_local, &pAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


		double alpha = rz_old / pAp;

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
		{
			omega[i] += alpha * p[i];

			r[i] -= alpha * Ap[i];

			z[i] = r[i] / D[i];
		}

		double rz_new_local = DotProduct(z, r);
		double rz_new = 0.0;
		MPI_Allreduce(&rz_new_local, &rz_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		double global_residual = std::sqrt(rz_new);
		if (world_rank == 0)
		{
			std::ofstream os("residual.txt",std::ios::app);
			os << global_residual << std::endl;
			os.close();
		}

		if (global_residual < delta)
		{
			if (world_rank == 0)
				std::cout << "converged in " << it << " stepts\n";

			break; // converged
		}

		double beta = rz_new / rz_old;

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
		{
			p[i] = z[i] + beta * p[i];
		}

		rz_old = rz_new;
	}

	if (world_rank == 0)
		std::cout << "steps limit (" << max_iter << ") reached" << std::endl;

	return omega;
}

void  MPINode::CreateDomainInfo(const Domain& InitialDomain, int M, int N)
{
	int world_rank, world_size;
	int X_segments, Y_segments;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	std::vector<Domain> Domains = domain::SplitDomain2D(world_size, InitialDomain);

	if (Domains.size() != world_size)
	{
		assert(false);
		MPI_Finalize();
		return;
	}

	m_Subdomain = Domains[world_rank];

	domain::FindOptimalPartitionRC(world_size, M, N, X_segments, Y_segments);

	BuildNeighborInfo(world_rank, X_segments, Y_segments);
}

void MPINode::BuildNeighborInfo(int world_rank, int X_segments, int Y_segments)
{
	m_vNeighbors.clear();
	m_vNeighbors.reserve(4);

	// Compute 2D grid coordinates of this rank
	int ix = world_rank % X_segments;
	int iy = world_rank / X_segments;
	int Nx_local = m_Subdomain.Nx_local;
	int Ny_local = m_Subdomain.Ny_local;
	int Nx_total = m_Subdomain.Nx_total;
	int Ny_total = m_Subdomain.Ny_total;

	auto addNeighbor = [&](int nx, int ny, NeighborDirection dir)
	{
		// Out of bounds = no neighbor
		if (nx < 0 || nx >= X_segments || ny < 0 || ny >= Y_segments)
		{
			m_aHasNeighbour[(int)dir] = false;
			return;
		}

		NeighborInfo nb;
		nb.NeighborRank = ny * X_segments + nx;
		nb.direction = dir;
		
		switch (dir)
		{
		case NeighborDirection::LEFT:
		{
			nb.sendIdx.reserve(Ny_local);
			nb.recvIdx.reserve(Ny_local);

			int y_begin = (iy == Y_segments - 1) ? 0 : 1;
			int y_end = y_begin + Ny_local;
			for (int y = y_begin; y < y_end; y++)
			{
				nb.sendIdx.push_back(y * Nx_total + 1);
				nb.recvIdx.push_back(y * Nx_total);
			}
		}
			break;
		case NeighborDirection::RIGHT:
		{
			nb.sendIdx.reserve(Ny_local);
			nb.recvIdx.reserve(Ny_local);

			int y_begin = (iy == Y_segments - 1) ? 0 : 1;
			int y_end = y_begin + Ny_local;
			for (int y = y_begin; y < y_end; y++)
			{
				nb.sendIdx.push_back((y + 1) * Nx_total - 2);
				nb.recvIdx.push_back((y + 1) * Nx_total - 1);
			}
		}
			break;
		case NeighborDirection::UP:
		{
			nb.sendIdx.reserve(Nx_local);
			nb.recvIdx.reserve(Nx_local);

			int x_begin = ix == 0 ? 0 : 1;
			int x_end = x_begin + Nx_local;
			for (int x = x_begin; x < x_end; x++)
			{
				nb.sendIdx.push_back(x + Nx_total);
				nb.recvIdx.push_back(x);
			}
		}
			break;
		case NeighborDirection::DOWN:
		{
			nb.sendIdx.reserve(Nx_local);
			nb.recvIdx.reserve(Nx_local);

			int x_begin = ix == 0 ? 0 : 1;
			int x_end = x_begin + Nx_local;
			for (int x = x_begin; x < x_end; x++)
			{
				nb.sendIdx.push_back(Nx_total * (Ny_total - 2) + x);
				nb.recvIdx.push_back(Nx_total * (Ny_total - 1) + x);
			}
		}
			break;
		}

		m_aHasNeighbour[(int)dir] = true;
		m_vNeighbors.push_back(nb);
	};

	// Add 4 possible neighbors
	addNeighbor(ix - 1, iy, NeighborDirection::LEFT);
	addNeighbor(ix + 1, iy, NeighborDirection::RIGHT);
	addNeighbor(ix, iy + 1, NeighborDirection::UP);
	addNeighbor(ix, iy - 1, NeighborDirection::DOWN);
}

int MPINode::GetHorizontalNeighboursCount()
{
	int count = 0;
	count += m_aHasNeighbour[(int)NeighborDirection::LEFT] ? 1 : 0;
	count += m_aHasNeighbour[(int)NeighborDirection::RIGHT] ? 1 : 0;
	return count;
}

int MPINode::GetVertialNeighboursCount()
{
	int count = 0;
	count += m_aHasNeighbour[(int)NeighborDirection::UP] ? 1 : 0;
	count += m_aHasNeighbour[(int)NeighborDirection::DOWN] ? 1 : 0;
	return count;
}

int MPINode::GetNeighboursCount()
{
	return GetHorizontalNeighboursCount() + GetVertialNeighboursCount();
}

void MPINode::GatherOmega(const std::vector<double>& omega, int M, int N, int X_segments, int Y_segments, bool bSave)
{
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Rank of the process
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int local_w = M + 1;
	int local_h = N + 1;
	int local_size = local_w * local_h;

	// --- Gather all sizes to root ---
	std::vector<int> recv_counts(world_size);
	MPI_Gather(&local_size, 1, MPI_INT,
		recv_counts.data(), 1, MPI_INT,
		0, MPI_COMM_WORLD);

	std::vector<int> displs(world_size);
	int offset = 0;

	if (world_rank == 0)
	{
		for (int i = 0; i < world_size; i++)
		{
			displs[i] = offset;
			offset += recv_counts[i];
		}
	}

	std::vector<double> gathered;
	if (world_rank == 0)
		gathered.resize(offset);

	// Gather ω from all ranks
	MPI_Gatherv(omega.data(), local_size, MPI_DOUBLE,
		gathered.data(), recv_counts.data(),
		displs.data(), MPI_DOUBLE,
		0, MPI_COMM_WORLD);


	if (world_rank == 0)
	{
		// Global grid resolution
		int global_w = X_segments * (local_w - 1) + 1;
		int global_h = Y_segments * (local_h - 1) + 1;

		std::vector<double> global_omega(global_w * global_h);

		// Reconstruct full grid
		for (int rank = 0; rank < world_size; rank++)
		{
			int i = rank % X_segments;   // column
			int j = rank / X_segments;   // row

			int gx0 = i * (local_w - 1);
			int gy0 = j * (local_h - 1);

			const double* block = &gathered[displs[rank]];

			int w_copy = (i == X_segments - 1) ? local_w : (local_w - 1);
			int h_copy = (j == Y_segments - 1) ? local_h : (local_h - 1);

			for (int y = 0; y < h_copy; y++)
			{
				for (int x = 0; x < w_copy; x++)
				{
					int gx = gx0 + x;
					int gy = gy0 + y;

					global_omega[gy * global_w + gx] = block[y * local_w + x];
				}
			}
		}

		if (bSave)
		{
			std::string ResultFileName = "Result" + std::to_string(M * X_segments) + "x" + std::to_string(N * Y_segments) + ".txt";
			std::ofstream ResultFile(ResultFileName);
			PrintFlatMatrix(ResultFile, global_omega, N * Y_segments + 1, M * X_segments + 1);
			ResultFile.close();
		}
	}
}