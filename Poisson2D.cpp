#include <iostream>
#include <limits>
#include <ctime>

#include "CSRMatrix.h"

constexpr double X_max = 1.0;
constexpr double Y_max = 1.0;
constexpr double X_min = -X_max;
constexpr double Y_min = -Y_max;


/*

		domain D
		   ↑ Y
 **********|----------
 **********|----------
 **********|----------
 **********|----------
 **********|---------- X
-----------|-----------→
 **********|**********
 **********|**********
 **********|**********
 **********|**********
 **********|**********
		   |

where '*' - inside D, '-' - outside

*/

inline bool IsEdge(double X, double Y)
{
	bool result = (std::abs(X) < FLT_EPSILON) && (Y > FLT_EPSILON);
	result = result || (std::abs(Y) < FLT_EPSILON) && (X > FLT_EPSILON);
	return result;
}

inline bool IsOutOfDomain(double X, double Y)
{
	return (X > FLT_EPSILON) && (Y > FLT_EPSILON);
}

void CreateMatrixes(CSRMatrix& A, std::vector<double>& f, __int64 M, __int64 N)
{
	__int64 Nn = N + 1, Mn = M + 1; // grid nodes count

	double X_step, Y_step, a = 0.0, b = 0.0;
	X_step = (X_max - X_min) / (double)N;
	Y_step = (Y_max - Y_min) / (double)M;
	double OneBy_h1 = 1.0 / (X_step * X_step); // 1.0/h_1^2 actually
	double OneBy_h2 = 1.0 / (Y_step * Y_step); // 1.0/h_2^2 actually
	double EPS = std::max(X_step * X_step, Y_step * Y_step);

	std::vector<Triplet> COO; // coordinate list matrix format
	COO.reserve(8 * Nn * Mn); // expected COO size

	f.clear();
	f.resize(Nn * Mn, 0.0);

	for (int i = N - 1; i > 0; i--) // from positive Y to negative
	{
		double Y_cur = Y_min + i * Y_step;
		for (int j = 1; j < M; j++) // from negative X to positive
		{
			__int64 idx = N - i;
			double X_cur = X_min + j * X_step;
			// default assumption - inside D
			a = b = 1.0;
			f[idx * Mn + j] = 1.0;

			if (IsEdge(X_cur, Y_cur)) // Edge
			{
				if (std::abs(X_cur) < EPS)
				{
					a = 0.5 * (1.0 / EPS + 1.0);
					b = 1.0;
				}
				else if (std::abs(Y_cur) < EPS)
				{
					a = 1.0;
					b = 0.5 * (1.0 / EPS + 1.0);
				}
				f[idx * Mn + j] = 0.5;
			}
			else if (IsOutOfDomain(X_cur, Y_cur)) // Outside
			{
				a = b = 1.0 / EPS;
				f[idx * Mn + j] = 0.0;
			}

			if (idx > 1)
			{
				COO.push_back({ (idx - 1) * Mn + j ,(idx - 1) * Mn + j,  a * OneBy_h1 });
				COO.push_back({ (idx - 1) * Mn + j ,idx * Mn + j,  -a * OneBy_h1 });
				COO.push_back({ idx * Mn + j ,(idx - 1) * Mn + j,  -a * OneBy_h1 });
			}

			if (j > 1)
			{
				COO.push_back({ idx * Mn + j - 1 ,idx * Mn + j - 1,  b * OneBy_h2 });
				COO.push_back({ idx * Mn + j - 1 ,idx * Mn + j,  -b * OneBy_h2 });
				COO.push_back({ idx * Mn + j ,idx * Mn + j - 1,  -b * OneBy_h2 });
			}

			COO.push_back({ idx * Mn + j ,idx * Mn + j,  a * OneBy_h1 + b * OneBy_h2 });
		}
	}

	f[Mn / 2 * Mn + Nn / 2] = 0.75; // center point 

	A = CSRMatrix::COO_To_CSR(COO, Nn * Mn, Nn * Mn);
}

std::vector<double> ConjugateGradient(const CSRMatrix& A, const std::vector<double>& F)
{
	int n = A.m_iRows;
	const int max_iter = n;
	const double delta = 1e-10;

	std::vector<double> omega(n, 0.0);
	std::vector<double> r = F;           // r0 = F - A*x = F
	std::vector<double> p = r;
	std::vector<double> Ap(n, 0.0);

	double rsold = 0.0;
	for (double v : r) 
		rsold += v * v;

	for (int it = 0; it < max_iter; it++) 
	{
		Ap = A.VectorMultiply(p);

		double pAp = 0.0;
		for (int i = 0; i < n; i++) 
			pAp += p[i] * Ap[i];

		double alpha = rsold / pAp;

		for (int i = 0; i < n; i++) 
		{
			omega[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		double rsnew = 0.0;
		for (double v : r) 
			rsnew += v * v;

		if (std::sqrt(rsnew) < delta) 
		{
			//std::cout << "Converged in " << it + 1 << " iterations\n";
			// converged
			break;
		}

		double beta = rsnew / rsold;
		for (int i = 0; i < n; i++) 
		{
			p[i] = r[i] + beta * p[i];
		}

		rsold = rsnew;
	}

	return omega;
}

int main()
{
	__int64 N, M; // X axis partitioned to M segments, Y - to N
	M = 40;
	N = 40;

	CSRMatrix A;
	std::vector<double> F, omega; // these are matrixes, just flatten
	CreateMatrixes(A, F, M, N);

	omega = ConjugateGradient(A, F);
}
