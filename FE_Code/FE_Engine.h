#ifndef FE_ENGINE_H
#define FE_ENGINE_H

#include <stdio.h>
#include <time.h>

#include "Array.h"
#include "Array2D.h"

#include "MatrixLibrary.h"

#include "FE_Element.h"
#include "FE_Mesh.h"
#include "Material.h"
#include "MatrixSolver.h"

#include "table_delaunay.h"

/**
 *  @class  FE_Engine
 *  @brief  Base class for a
 */
class FE_Engine {

public:
	FE_Engine();

	FE_Engine(FILE* File);

	FE_Engine(int nx,
			  int ny,
			  double xmin,
			  double xmax,
			  double ymin,
			  double ymax,
			  double xscale,
			  double yscale);

	virtual ~FE_Engine();

	void compute_stiffness_matrix(	Material &Material_,
									FE_ElementQuad &feElement_);

	void compute_rhs_vector(DENS_VEC &bc,
							DENS_VEC &value);

	void compute_dirichlet_bc(DENS_VEC &bc,
								DENS_VEC &value,
								DENS_VEC &F);

	void compute_dirichlet_bc_y(DENS_VEC &bc,
								DENS_VEC &value,
								DENS_VEC &F);

	void compute_dirichlet_bc_x(DENS_VEC &bc,
								DENS_VEC &value,
								DENS_VEC &F);

	void compute_neumann_bc() const;

	void driver(double mu,
				double nu,
				int ndof,
				int ncoord,
				double load,
				DENS_VEC dirichletprescnodes,
				DENS_VEC dirichletpresc_values,
				DENS_VEC &updatedCoords,
				DENS_VEC &lambda) const;

	void twod_to_vtk (int tsp, DENS_MAT &xyz,
					  Array2D<int> &e_conn,
					  DENS_VEC &p,
					  DENS_MAT &uvw);

	void vtk_puv_write (FILE* output_unit,
						int nNodes,
						int nElements,
						int element_order,
						DENS_MAT &xyz,
						Array2D<int> &e_conn,
						DENS_VEC &p, DENS_MAT &uvw ) const;

	double diffclock(clock_t clock1,clock_t clock2) const;

	void post_processing(int tsp, DENS_VEC &lambda);

	void compute_dirichlet_presc(DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F, DENS_MAT &K);
	
	void compute_robin(DENS_MAT Neighbor, DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F, DENS_MAT &K, double rcoeff);

	// before Robin
	// void FESolve(int tsp, DENS_VEC drchlt_prsc_nds,
	//		DENS_VEC drchlt_prsc_vls,
	//		DENS_VEC &updatedCoords,
	//		DENS_VEC &lambda, int flag);
	
	
	void FESolve(int tsp, DENS_VEC drchlt_prsc_nds,
			DENS_VEC drchlt_prsc_vls,
			DENS_VEC &updatedCoords,
			DENS_VEC &lambda, int flag, 
            double rcoeff, DENS_MAT Neighbor);

	void FEsetup(double mu,
			double nu,
			int ndof,
			int ncoord,
			double load);
	
	void FEsetup1(double mu,
			double nu,
			int ndof,
			int ncoord,
			double load);

	void FEsetupload1(double mu,
			double nu,
			int ndof,
			int ncoord,
			double load);

	// used for Neumann (June 27/2016)
	void GenerateMesh(DENS_MAT coords_xy0, DENS_MAT &Connectivity);
	void ComputeStress(double PoissonRatio, double YoungModulus, 
				      DENS_MAT coords_xy0, DENS_MAT coords_xy,
				      DENS_MAT Connectivity, DENS_MAT &Stress);

	template <typename T>
	std::string number_to_string(T number);

	// Mesh obj
	FE_Uniform2DMesh * feMesh;

protected:

//	// Mass Matrix
//	SPAR_MAT mass_;
//
//	// Stiffness Matrix
//	SPAR_MAT stiff_;
//
//	// Right hand side
//	DENS_MAT rhs_;

	int nx_;
	int ny_;
	double xmin_;
	double xmax_;
	double ymin_;
	double ymax_;
	double xscale_;
	double yscale_;

	// Stiffness Matrix
	DENS_MAT GlbStiffMatrix;

	// Right hand side
	DENS_VEC rhs;

};

#endif // FE_ENGINE_H
