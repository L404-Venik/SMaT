#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <omp.h>


struct PairHash 
{
	std::size_t operator()(const std::pair<int, int>& p) const noexcept 
	{
		return (static_cast<std::size_t>(p.first) << 32) ^ static_cast<std::size_t>(p.second);
	}
};

struct Triplet
{
	int row, col;
	double value;
};

class CSRMatrix // Compressed sparse row Matrix class
{
public:
	int m_iRows, m_iCols;
	std::vector<double> values;
	std::vector<int> col_index;
	std::vector<int> row_ptr;

	CSRMatrix() :m_iRows(0), m_iCols(0) {};
	CSRMatrix(int r, int c);

	// Convert COO → CSR
	static CSRMatrix COO_To_CSR(const std::vector<Triplet>& coo, int rows, int cols);
	
	// Sparse matrix-vector multiply: y = A * x
	std::vector<double> VectorMultiply(const std::vector<double>& x) const;
	std::vector<double> GetDiagonal() const;

	bool is_symmetric(double tol = 1e-12) const;
	void print(std::ofstream& OutStream) const;
};


/**
* @brief Computes dot product of 2 vectors of the same size
* @return Dot product value.
*/
double DotProduct(const std::vector<double>& x, const std::vector<double>& y);

/**
* @brief Prints flattened matrix to a file
*
* @param sFilePath - Path to a file to print matrix into.
* @param matrix - Flattened matrix.
* @param Nx - Rows number.
* @param Ny - Columns number.
*/
void PrintFlatMatrix(const std::string& sFilePath, const std::vector<double>& matrix, int Nx, int Ny);