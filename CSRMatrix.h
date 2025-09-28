#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

struct Triplet 
{
    __int64 row, col;
    double value;
};

class CSRMatrix // Compressed sparse row Matrix class
{
public:
    __int64 m_iRows, m_iCols;
    std::vector<double> values;
    std::vector<__int64> col_index;
    std::vector<__int64> row_ptr;

    CSRMatrix() :m_iRows(0), m_iCols(0) {};
    CSRMatrix(__int64 r, __int64 c);

    // Convert COO → CSR
    static CSRMatrix COO_To_CSR(const std::vector<Triplet>& coo, int rows, int cols);

    // Sparse matrix-vector multiply: y = A * x
    std::vector<double> VectorMultiply(const std::vector<double>& x) const;
};