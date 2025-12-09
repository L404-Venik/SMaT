#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <cstring>

#include "CSRMatrix.h"

int triplet_cmp(const void* a, const void* b)
{
	const Triplet* ta = static_cast<const Triplet*>(a);
	const Triplet* tb = static_cast<const Triplet*>(b);

	if (ta->row < tb->row) return -1;
	if (ta->row > tb->row) return 1;
	if (ta->col < tb->col) return -1;
	if (ta->col > tb->col) return 1;
	return 0;
}

// Constructors
CSRMatrix::CSRMatrix(int r, int c) : m_iRows(r), m_iCols(c)
{
	row_ptr.resize(r + 1, 0);

	// this Poisson task specific size requirements
	col_index.reserve(5 * r);
	values.reserve(5 * r);
}


// Convert COO → CSR
CSRMatrix CSRMatrix::COO_To_CSR(const std::vector<Triplet>& coo, int rows, int cols)
{
	CSRMatrix A(rows, cols);

	// Copy and sort by (row, col)
	std::vector<Triplet> sorted = coo;
	qsort(sorted.data(), sorted.size(), sizeof(Triplet), triplet_cmp);

	// Merge duplicates and fill CSR
	int nnz = 0;
	A.row_ptr[0] = 0;
	int current_row = 0;

	for (size_t k = 0; k < sorted.size(); )
	{
		int i = sorted[k].row;
		int j = sorted[k].col;
		double val = 0.0;

		// accumulate duplicates (same row & col)
		while (k < sorted.size() && sorted[k].row == i && sorted[k].col == j)
		{
			val += sorted[k].value;
			k++;
		}

		// fill skipped rows with row_ptr update
		while (current_row < i)
		{
			A.row_ptr[current_row + 1] = nnz;
			current_row++;
		}

		// store element
		A.values.push_back(val);
		A.col_index.push_back(j);
		nnz++;
	}

	// finalize row_ptr
	while (current_row < rows)
	{
		A.row_ptr[current_row + 1] = nnz;
		current_row++;
	}

	return A;
}

// Sparse matrix-vector multiply: y = A * x
std::vector<double> CSRMatrix::VectorMultiply(const std::vector<double>& x) const
{
	if (x.size() != m_iCols)
		throw std::runtime_error("Dimension mismatch in mat-vec multiplication");

	std::vector<double> y(m_iRows, 0.0);

	#pragma omp parallel for schedule(static)// num_threads(NumThreads)
	for (int i = 0; i < m_iRows; i++)
	{
		double sum = 0.0;
		int start = row_ptr[i];
		int end = row_ptr[i + 1];
		for (int k = start; k < end; k++)
			sum += values[k] * x[col_index[k]];

		y[i] = sum;
	}

	return y;
}

std::vector<double> CSRMatrix::GetDiagonal() const
{
	std::vector<double> D(m_iRows, 1.0);

	for (int i = 0; i < m_iRows; ++i)
	{
		for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
		{
			if (col_index[k] == i)
			{
				D[i] = values[k];
			}
		}
	}

	return D;
}

// Other funtions

bool CSRMatrix::is_symmetric(double tol) const
{
	if (m_iRows != m_iCols)
		return false;

	for (int i = 0; i < m_iRows; i++)
	{
		for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
		{
			int j = col_index[k];
			double a_ij = values[k];

			// Find A[j,i]
			bool found = false;
			for (int kk = row_ptr[j]; kk < row_ptr[j + 1]; kk++)
			{
				if (col_index[kk] == i)
				{
					if (std::fabs(values[kk] - a_ij) > tol) return false;
					found = true;
					break;
				}
			}
			if (!found && std::fabs(a_ij) > tol)
				return false;
		}
	}

	return true;
}

void CSRMatrix::print(std::ofstream& OutStream) const
{
	for (int i = 0; i < m_iRows; i++)
	{
		int start = row_ptr[i];
		int end = row_ptr[i + 1];

		// fill row with zeros
		std::vector<double> row(m_iCols, 0.0);

		// place nonzeros
		for (int k = start; k < end; k++)
		{
			row[col_index[k]] = values[k];
		}

		// print row
		for (int j = 0; j < m_iCols; j++)
		{
			OutStream << std::setw(8) << row[j] << " ";
		}
		OutStream << "\n";
	}
}

//////////////////////////////////////////////////////////////////////////////
// Other common functions

double DotProduct(const std::vector<double>& x, const std::vector<double>& y)
{
	assert(x.size() == y.size());

	double result = 0.0;
	size_t n = x.size();

	#pragma omp parallel for reduction(+:result) schedule(static)// num_threads(NumThreads)
	for (int i = 0; i < n; i++)
		result += x[i] * y[i];

	return result;
}

// Nx - x nodes count, Ny - y nodes count
void PrintFlatMatrix(std::ofstream& output, const std::vector<double>& matrix, int Nx, int Ny)
{
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			output << matrix[i * Ny + j] << ' ';
		}
		output << std::endl;
	}
}

void PrintFlatMatrix(const std::string& sFilePath, const std::vector<double>& matrix, int Nx, int Ny)
{
	std::ofstream output(sFilePath);
	PrintFlatMatrix(output, matrix, Nx, Ny);
	output.close();
}