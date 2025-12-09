#pragma once
#include <algorithm>
#include "CSRMatrix.h"

/**
 * @brief Solves linear equations system with conjugate gradien method.
 *
 * @param A - CSRMatrxi of linear algebraic equations equations.
 * @param F - Right-hand side vector.
 * @param delta - Euclidean norm of discrepancy threshold.
 * @return std::vector<double> - approximate system solution.
 */
std::vector<double> ConjugateGradient(const CSRMatrix& A, const std::vector<double>& F, double delta = 0.05);