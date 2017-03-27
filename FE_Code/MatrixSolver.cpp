#include <iostream>
using namespace std;

#include "MatrixSolver.h"
#include "CG.h"
//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolver
//--------------------------------------------------------
//--------------------------------------------------------

//--------------------------------------------------------
//  Constructor
//         Grab references to necessary data
//--------------------------------------------------------
LambdaMatrixSolver::LambdaMatrixSolver(int maxIterations, double tolerance) :
maxIterations_(maxIterations),
tolerance_(tolerance)
{
	// do nothing
}

LambdaMatrixSolver::~LambdaMatrixSolver()
{

}

//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolverLumped
//--------------------------------------------------------
//--------------------------------------------------------

//--------------------------------------------------------
//  Constructor
//         Grab references to necessary data
//--------------------------------------------------------
LambdaMatrixSolverLumped::LambdaMatrixSolverLumped(int maxIterations, double tolerance) :
    		LambdaMatrixSolver(maxIterations,tolerance)
{
	// do nothing
}

void LambdaMatrixSolverLumped::execute(DENS_VEC & rhs, DENS_VEC & lambda, DIAG_MAT & lumpedMatrix)
{
	// solve lumped equation
	for (int i = 0; i < rhs.size(); i++)
		lambda(i) = rhs(i)/lumpedMatrix(i,i);
}

LambdaMatrixSolverLumped::~LambdaMatrixSolverLumped()
{

}
//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolverCg
//--------------------------------------------------------
//--------------------------------------------------------

//--------------------------------------------------------
//  Constructor
//         Grab references to necessary data
//--------------------------------------------------------
LambdaMatrixSolverCg::LambdaMatrixSolverCg(int maxIterations, double tolerance) :
    		LambdaMatrixSolver(maxIterations, tolerance)
{
	// do nothing
}

void LambdaMatrixSolverCg::execute(DENS_VEC & rhs, DENS_VEC & lambda, DENS_MAT & myMatrix)
{
	DIAG_MAT preConditioner = myMatrix.get_diag();
	int myMaxIt  = 2*myMatrix.nRows(); // note could also use the fixed parameter
	double myTol = tolerance_;

	cout << " " << endl;
	cout << " myMaxIt...:" << " " << myMaxIt << endl;
	cout << " myTol.....:" << " " << myTol << endl;

	int convergence = CG(myMatrix, lambda, rhs, preConditioner, myMaxIt, myTol);

	cout << " " << endl;
	cout << " Converged at...:" << " " << myMaxIt << endl;
	cout << " Residual.......:" << " " << myTol << endl;

	// error if didn't converge
	if (convergence>0)
		cout << "CG solver did not converge in LambdaMatrixSolverCg::execute()" << endl;
}


LambdaMatrixSolverCg::~LambdaMatrixSolverCg()
{
	// do nothing
}
