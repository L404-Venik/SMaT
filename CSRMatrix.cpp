#include <omp.h>
#include "CSRMatrix.h"

CSRMatrix::CSRMatrix(__int64 r, __int64 c) : m_iRows(r), m_iCols(c)
{
    row_ptr.resize(r + 1, 0);

    // this Poisson task specific size requirements
    col_index.reserve(5 * r);
    values.reserve(5 * r);
}

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

    //#pragma omp parallel for
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