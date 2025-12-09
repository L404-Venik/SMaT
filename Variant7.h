#pragma once
#include "CSRMatrix.h"
#include "Domain.h"

/*  L-shaped Domain, axis rotated 90 degrees clockwise
	{(x, y) : −1 < x, y < 1} \ {(x, y) : 0 < x, y < 1}
		domain D  
		   -1.0|
	 **********|**********
	 **********|**********
	 **********|**********
	 **********|**********
	 **********|********** Y
   ------------0-----------→
-1.0 **********|---------- 1.0
	 **********|----------
	 **********|----------
	 **********|----------
	 **********|----------
			1.0↓ X

where '*' - inside D, '-' - outside

*/

constexpr double X_max = 1.0;
constexpr double Y_max = 1.0;
constexpr double X_min = -X_max;
constexpr double Y_min = -Y_max;

/**
* @brief Create matrix and right-hand side vector of linear equations system for 2D Poisson problem in L-shaped domain.
*
* @param A - Destination for a CSRMatrix reference.
* @param F - Destination for a right-hand side vector reference.
* @param Nx - Number of nodes by X axis.
* @param Ny - Number of nodes by Y axis.
* @param D - Domain to make matrixes for.
*/
void CreateMatrixesV7(CSRMatrix& A, std::vector<double>& F, int Nx, int Ny, const Domain& D = { X_min, X_max, Y_min, Y_max });