#include <iostream>
#include <limits>
#include <ctime>

constexpr double X_min = -1;
constexpr double Y_min = -1;

constexpr double X_max = 1;
constexpr double Y_max = 1;

/*

	domain D
**********----------
**********----------
**********----------
**********----------
**********----------
********************
********************
********************
********************
********************

where '*' - inside D, '-' - outside

*/

inline bool IsEdge(double X, double Y)
{
	bool result = (std::abs(X) < FLT_EPSILON) && (Y > FLT_EPSILON);
	result = result || (std::abs(Y) < FLT_EPSILON) && (X > FLT_EPSILON);
	return result;
}

inline bool IsOutOfBody(double X, double Y)
{
	return (X > FLT_EPSILON) && (Y > FLT_EPSILON);
}

int main()
{
	//double start = std::qua

	__int64 N, M; // X axis partitioned to M segments, Y - to N
	N = M = 10;

	double X_step, Y_step, a = 0.0, b = 0.0;
	X_step = (X_max - X_min) / (double)N;
	Y_step = (Y_max - Y_min) / (double)M;
	double OneBy_h1 = 1.0 / (X_step * X_step); // 1.0/h_1^2 actually
	double OneBy_h2 = 1.0 / (Y_step * Y_step); // 1.0/h_2^2 actually

	double* A_flat = new double[N * M * N * M];
	memset(A_flat, 0.0, N * M * N * M * sizeof(double));

	double** A = new double*[N * M];
	for (int i = 0; i < N * M; ++i)
		A[i] = &A_flat[i * N * M];

	double* d_omega = new double[N * M];
	memset(d_omega, 0.0, N * M * sizeof(double));

	double* f = new double[N * M];
	memset(f, 1.0, N * M * sizeof(double));

	for (int i = 1; i < N - 1; i++)
	{
		double X_cur = X_min + i * X_step;
		for (int j = 1; j < M - 1; j++)
		{
			double Y_cur = Y_min + i * Y_step;
			if (IsEdge(X_cur, Y_cur))
			{
				if (std::abs(X_cur) < FLT_EPSILON)
				{
					a = 0.5 * (1.0 / FLT_EPSILON + 1.0);
					b = 1.0;
				}
				else if (std::abs(Y_cur) < FLT_EPSILON)
				{
					a = 1.0;
					b = 0.5 * (1.0 / FLT_EPSILON + 1.0);
				}
				f[i * M + j] = 0.5;

			}
			else if (IsOutOfBody(X_cur, Y_cur))
			{
				a = b = 1.0 / FLT_EPSILON;
				f[i * M + j] = 0.0;
			}
			else
			{
				f[i * M + j] = 1.0;

				//if(std::abs(X_cur) < FLT_EPSILON && std::abs(Y_cur) < FLT_EPSILON) // point (0,0)
				//	f[i * M + j] = 0.75;

				a = b = 1.0;
			}

			A[i * M + j][i * M + j] += a * OneBy_h1;
			A[(i - 1) * M + j][(i - 1) * M + j] += a * OneBy_h1;
			A[(i - 1) * M + j][i * M + j] -= a * OneBy_h1;
			A[i * M + j][(i - 1) * M + j] -= a * OneBy_h1;

			A[i * M + j][i * M + j] += b * OneBy_h2;
			A[i * M + j - 1][i * M + j - 1] += b * OneBy_h2;
			A[i * M + j - 1][i * M + j] -= b * OneBy_h2;
			A[i * M + j][i * M + j - 1] -= b * OneBy_h2;
		}
	}

	f[ M / 2 * M + N / 2] = 0.75;

	delete[] d_omega;
	delete[] f;
	delete[] A_flat;
}
