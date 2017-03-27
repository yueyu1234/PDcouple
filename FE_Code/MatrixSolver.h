#ifndef MATRIXSOLVER_H
#define MATRIXSOLVER_H

#include "MatrixLibrary.h"

//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolver
//--------------------------------------------------------
//--------------------------------------------------------

class LambdaMatrixSolver {

public:

	LambdaMatrixSolver(int maxIterations, double tolerance);

	~LambdaMatrixSolver();

	/** execute the solver */
	//virtual void execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & lumpedMatrix) = 0;

protected:

	/** maximum number of iterations */
	int maxIterations_;

	/** relative tolerance to solve to */
	double tolerance_;

private:

	// DO NOT define this
	LambdaMatrixSolver();

};

//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolverLumped
//--------------------------------------------------------
//--------------------------------------------------------

class LambdaMatrixSolverLumped : public LambdaMatrixSolver {

public:

	LambdaMatrixSolverLumped(int maxIterations, double tolerance);

	~LambdaMatrixSolverLumped();

	/** execute the solver */
	void execute(DENS_VEC & rhs, DENS_VEC & lambda, DIAG_MAT & lumpedMatrix);

protected:


private:

	// DO NOT define this
	//LambdaMatrixSolverLumped();

};

//--------------------------------------------------------
//--------------------------------------------------------
//  Class LambdaMatrixSolverCg
//--------------------------------------------------------
//--------------------------------------------------------

class LambdaMatrixSolverCg : public LambdaMatrixSolver {

public:

	LambdaMatrixSolverCg(int maxIterations, double tolerance);

	~LambdaMatrixSolverCg();

	/** execute the solver */
	void execute(DENS_VEC & rhs, DENS_VEC & lambda, DENS_MAT & myMatrix);

protected:



private:

	// DO NOT define this
	//LambdaMatrixSolverCg();

};

#endif // MATRIXSOLVER_H
