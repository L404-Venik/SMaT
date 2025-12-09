#include "Variant7.h"
#include <numeric>
#include <limits>
#include <float.h>
#include <cmath>
#include <cassert>
#include <string>

inline bool IsCenter(double X, double Y)
{
	return (std::abs(X) < FLT_EPSILON) && (std::abs(Y) < FLT_EPSILON);
}

inline bool IsInnerCorner(double X, double Y)
{
	return (std::abs(X) < FLT_EPSILON) && (Y > FLT_EPSILON) || (std::abs(Y) < FLT_EPSILON) && (X > FLT_EPSILON);
}

inline bool IsTopPart(double X, double Y)
{
	return std::abs(Y - Y_max) < FLT_EPSILON && X <= 0.0;
}

inline bool IsBottomPart(double X, double Y)
{
	return std::abs(Y - Y_min) < FLT_EPSILON && (X <= 1.0 && X >= -1.0);
}

inline bool IsLeftPart(double X, double Y)
{
	return std::abs(X - X_min) < FLT_EPSILON && (Y <= 1.0 && Y >= -1.0);
}

inline bool IsRightPart(double X, double Y)
{
	return std::abs(X - X_max) < FLT_EPSILON && Y <= 0.0;
}

inline bool IsEdge(double X, double Y)
{
	bool result = IsInnerCorner(X, Y)
		|| IsTopPart(X, Y)
		|| IsRightPart(X, Y)
		|| IsLeftPart(X, Y)
		|| IsBottomPart(X, Y);
	return result;
}

inline bool IsOuterCorner(double X, double Y)
{
	bool result = false;
	result = result || std::abs(X - X_min) < FLT_EPSILON && std::abs(Y - Y_max) < FLT_EPSILON; // left top
	result = result || std::abs(X) < FLT_EPSILON && std::abs(Y - Y_max) < FLT_EPSILON; // top center
	result = result || std::abs(X - X_max) < FLT_EPSILON && std::abs(Y) < FLT_EPSILON; // center right
	result = result || std::abs(X - X_max) < FLT_EPSILON && std::abs(Y - Y_min) < FLT_EPSILON; // right bottom
	result = result || std::abs(X - X_min) < FLT_EPSILON && std::abs(Y - Y_min) < FLT_EPSILON; // bottom left
	return result;
}

inline bool IsOutOfDomain(double X, double Y)
{
	return (X > FLT_EPSILON) && (Y > FLT_EPSILON);
}

// Nx - x nodes count, Ny - y nodes count
void CreateWithCOO(CSRMatrix& A, std::vector<double>& F, int Nx, int Ny, const Domain& D)
{
	double X_step, Y_step, a = 0.0, b = 0.0;
	X_step = (D.x_max - D.x_min) / (double)(Nx - 1);
	Y_step = (D.y_max - D.y_min) / (double)(Ny - 1);
	double OneBy_h1 = 1.0 / (X_step * X_step); // 1.0/h_1^2 actually
	double OneBy_h2 = 1.0 / (Y_step * Y_step); // 1.0/h_2^2 actually
	double EPS = std::max(X_step * X_step, Y_step * Y_step);

	std::vector<double> A_mat(Nx * Ny, 0.0), B_mat(Nx * Ny, 0.0);
	std::vector<Triplet> COO; // coordinate list matrix format
	COO.reserve(8 * Ny * Nx); // expected COO size

	F.clear();
	F.resize(Ny * Nx, 0.0);

	for (int i = 0; i < Nx; i++) // from Negative X to positive
	{
		double X_cur = D.x_min + i * X_step;

		for (int j = 0; j < Ny; j++) // from negative Y to positive
		{
			double Y_cur = D.y_min + j * Y_step;
			// default assumption - inside D
			a = b = 1.0;
			F[i * Ny + j] = 1.0;

			if (IsEdge(X_cur, Y_cur)) // Edge
			{
				F[i * Ny + j] = 0.5;

				if (IsOuterCorner(X_cur, Y_cur))
				{
					a = b = 0.5 * (1.0 / EPS + 1.0);
					F[i * Ny + j] = 0.25;
				}
				else if (IsTopPart(X_cur, Y_cur) || IsBottomPart(X_cur, Y_cur) || (std::abs(Y_cur) < EPS && !IsLeftPart(X_cur, Y_cur))) // horizontal parts of edge
				{
					a = 0.5 * (1.0 / EPS + 1.0);
					b = 1.0;
				}
				else if (IsLeftPart(X_cur, Y_cur) || IsRightPart(X_cur, Y_cur) || (std::abs(X_cur) < EPS && !IsBottomPart(X_cur, Y_cur))) // vertical parts of edge
				{
					a = 1.0;
					b = 0.5 * (1.0 / EPS + 1.0);
				}

				if (IsLeftPart(X_cur, Y_cur))
				{
					a = 1 / EPS;
				}
				if (IsBottomPart(X_cur, Y_cur))
				{
					b = 1 / EPS;
				}

			}
			else if (IsOutOfDomain(X_cur, Y_cur)) // Outside
			{
				a = b = 1.0 / EPS;
				F[i * Ny + j] = 0.0;
			}
			else if (IsCenter(X_cur, Y_cur))
			{
				F[i * Ny + j] = 0.75;
			}

			if (i > 0)
			{
				COO.push_back({ (i - 1) * Ny + j ,(i - 1) * Ny + j,  a * OneBy_h1 });
				COO.push_back({ (i - 1) * Ny + j ,i * Ny + j,  -a * OneBy_h1 });
				COO.push_back({ i * Ny + j ,(i - 1) * Ny + j,  -a * OneBy_h1 });
			}

			if (j > 0)
			{
				COO.push_back({ i * Ny + j - 1 ,i * Ny + j - 1,  b * OneBy_h2 });
				COO.push_back({ i * Ny + j - 1 ,i * Ny + j,  -b * OneBy_h2 });
				COO.push_back({ i * Ny + j ,i * Ny + j - 1,  -b * OneBy_h2 });
			}

			COO.push_back({ i * Ny + j ,i * Ny + j,  a * OneBy_h1 + b * OneBy_h2 });

			A_mat[i * Ny + j] = a;
			B_mat[i * Ny + j] = b;
		}
	}

	/*std::string AFileName = "Amat.txt";
	PrintFlatMatrix(AFileName, A_mat, Ny, Nx);
	std::string BFileName = "Bmat.txt";
	PrintFlatMatrix(BFileName, B_mat, Ny, Nx);
	std::string FFileName = "Fmat.txt";
	PrintFlatMatrix(FFileName, F, Ny, Nx);*/

	A = CSRMatrix::COO_To_CSR(COO, Ny * Nx, Ny * Nx);
}

void CreateMatrixesV7(CSRMatrix& A, std::vector<double>& F, int Nx, int Ny, const Domain& D) // Variant 7
{
	CreateWithCOO(A, F, Nx, Ny, D);
}