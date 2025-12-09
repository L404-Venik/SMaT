#include <cstring>

#include "ConjugateGradient.h"

std::vector<double> CustomRealization(const CSRMatrix& A, const std::vector<double>& F, double delta)
{
	int n = A.m_iRows;
	const int max_iter = n;

	std::vector<double> omega(n, 0.0);
	std::vector<double> r = F;          // r0 = F - A*x = F
	std::vector<double> p;
	std::vector<double> Ap(n, 0.0), z(n, 0.0);
	std::vector<double> D = A.GetDiagonal();

	std::ofstream os("residual.txt", std::ios::out);// clear residual file
	os.close();

	for (int i = 0; i < n; ++i)
		z[i] = r[i] / D[i];

	p = z;
	double rz_old = DotProduct(z, r);

	for (int it = 0; it < max_iter; it++)
	{
		Ap = A.VectorMultiply(p);

		double pAp = DotProduct(p, Ap);

		double alpha = rz_old / pAp;

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
		{
			omega[i] += alpha * p[i];

			r[i] -= alpha * Ap[i];

			z[i] = r[i] / D[i];
		}

		double rz_new = DotProduct(z, r);
		double residual = std::sqrt(rz_new);

		std::ofstream os("residual.txt", std::ios::app);
		os << residual << std::endl;
		os.close();

		if (residual < delta)
		{
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

	return omega;
}

std::vector<double> StdRealization(const CSRMatrix& A, const std::vector<double>& F, double delta)
{
	int n = A.m_iRows;
	int M = std::sqrt(n);
	const int max_iter = n;

	std::vector<double> omega(n, 0.0);
	std::vector<double> r = F;          // r0 = F - A*x = F
	std::vector<double> p;
	std::vector<double> Ap(n, 0.0), z(n, 0.0);
	std::vector<double> D = A.GetDiagonal();
	std::vector<double> r_norms;

	for (int i = 0; i < n; ++i)
		z[i] = r[i] / D[i];

	p = z;
	double rz_old = DotProduct(z, r);

	for (int it = 0; it < max_iter; it++)
	{
		Ap = A.VectorMultiply(p);

		double pAp = std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);

		double alpha = rz_old / pAp;


		std::transform(p.begin(), p.end(), omega.begin(), omega.begin(),
			[alpha](double pi, double oi) { return oi + alpha * pi; });

		std::transform(Ap.begin(), Ap.end(), r.begin(), r.begin(),
			[alpha](double Ap, double r) { return r - alpha * Ap; });

		std::transform(r.begin(), r.end(), D.begin(), z.begin(),
			[](double r, double D) { return r / D; });

		double rz_new = std::inner_product(z.begin(), z.end(), r.begin(), 0.0);

		if (std::sqrt(rz_new) < delta)
		{
			std::cout << "converged in " << it << " stepts\n";
			break; // converged
		}
		//r_norms.push_back(rz_new);

		double beta = rz_new / rz_old;
		std::transform(z.begin(), z.end(), p.begin(), p.begin(),
			[beta](double zi, double pi) { return zi + beta * pi; });

		rz_old = rz_new;
	}

	/*std::ofstream ResultFile("residual.txt");
	PrintFlatMatrix(ResultFile, r_norms, r_norms.size(),1);
	ResultFile.close();*/

	return omega;
}

std::vector<double> ConjugateGradient(const CSRMatrix& A, const std::vector<double>& F, double delta)
{
	return CustomRealization(A, F, delta);
	//return StdRealization(A, F);
}