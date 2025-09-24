#include <iostream>
#include <limits>

constexpr double X_min = -1;
constexpr double Y_min = -1;

constexpr double X_max = 1;
constexpr double Y_max = 1;

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
	int N, M; // X axis partitioned to M segments, Y - to N
	N = M = 10;

	double X_step, Y_step, a, b;
	X_step = (X_max - X_min) / (double)N;
	Y_step = (Y_max - Y_min) / (double)M;

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
					a = 0.5 * (1.0 / FLT_EPSILON + 1);
					b = 1.0;
				}
				else if (std::abs(Y_cur) < FLT_EPSILON)
				{
					a = 1.0;
					b = 0.5 * (1.0 / FLT_EPSILON + 1);
				}
			}
			else if (IsOutOfBody(X_cur, Y_cur))
			{
				a = b = 1.0 / FLT_EPSILON;
				f[i * M + j] = 0;
			}
			else
			{
				a = b = 1.0;
			}
		}
	}


	std::pair<double, double>* Mesh = new std::pair<double, double>[N * M];



	delete[] Mesh;
	delete[] d_omega;
	delete[] f;
}
