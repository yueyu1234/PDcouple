#include "FE_Engine.h"

#include <iostream>
#include <cstdlib>
using namespace std;

#include <string>
#include <sstream>

// Other headers
#include "math.h"


// -------------------------------------------------------------
// -------------------------------------------------------------
//   class FE_Engine
// -------------------------------------------------------------
// -------------------------------------------------------------
FE_Engine::FE_Engine()
{
	feMesh = NULL;
}

FE_Engine::FE_Engine(FILE* File)
{
	feMesh = new FE_Uniform2DMesh(File);
}

FE_Engine::FE_Engine( int nx,
		int ny,
		double xmin,
		double xmax,
		double ymin,
		double ymax,
		double xscale,
		double yscale)
:nx_(nx),
 ny_(ny),
 xmin_(xmin),
 xmax_(xmax),
 ymin_(ymin),
 ymax_(ymax),
 xscale_(xscale),
 yscale_(yscale)
{
	//	fprintf ( stdout, "\n Generating mesh... " );
	feMesh = new FE_Uniform2DMesh(nx_, ny_, xmin_, xmax_, ymin_, ymax_, xscale_, yscale_);
}

FE_Engine::~FE_Engine()
{
	//if (feMesh)  delete feMesh;
}


//-----------------------------------------------------------------
void FE_Engine::compute_stiffness_matrix(Material &Material_, FE_ElementQuad &feElement_)
{
	int rw, cl;
	int nElements = feMesh->get_nElements();
	int nNodes    = feMesh->get_nNodes();

	int nNodesPerElement = feElement_.num_elt_nodes();
	int nDofsPerNode = feElement_.num_node_dofs();
	int nTotalDofs = nDofsPerNode * nNodes;

	//    fprintf ( stdout, "\n Calculation material matrix... " );
	// 1: PlaneStress
	double **** C = Material_.material_stiffness(1);

	// assemble consistent mass (nnodes X nnodes)
	GlbStiffMatrix.reset(nTotalDofs, nTotalDofs);
	DENS_MAT elemstiff(nDofsPerNode * nNodesPerElement, nDofsPerNode * nNodesPerElement);
	Array<int> conn(nNodesPerElement);
	DENS_MAT coords;
	for(int ielem = 0; ielem < nElements; ++ielem) {

		feMesh->element_coordinates(ielem, coords);
		feElement_.element_stiffness(C, coords, elemstiff);

		//		cout << " " << endl;
		//		elemstiff.print();
		//		cout << " " << endl;

		// get connectivity
		feMesh->element_connectivity_global(ielem, conn);

		// assembly
		for (int a = 0; a < nNodesPerElement; ++a){
			int inode = conn(a);
			for(int i = 0; i < nDofsPerNode; ++i) {
				for (int b = 0; b < nNodesPerElement; ++b) {
					int jnode = conn(b);
					for(int k = 0; k < nDofsPerNode; ++k) {
						rw = nDofsPerNode * inode + i;
						cl = nDofsPerNode * jnode + k;
						GlbStiffMatrix( rw, cl ) = ( GlbStiffMatrix( rw, cl ) + elemstiff(nDofsPerNode*a+i, nDofsPerNode*b+k) );
						//						matrix.add( rw, cl, elemstiff(nDofsPerNode*a+i, nDofsPerNode*b+k) );
					}
				}
			}
		}
		elemstiff.zero();
	}
	//		for (int i = 0; i < nNodesPerElement; ++i)
	//		{
	//			int inode = conn(i);
	//			for (int j = 0; j < nNodesPerElement; ++j)
	//			{
	//				int jnode = conn(j);
	//				matrix.add(inode, jnode, elemstiff(i,j));
	//			}
	//		}

	//	for a = 1:nelnodes(lmn)
	//	   for i = 1:nDofsPerNode
	//	       for b = 1:nelnodes(lmn)
	//	          for k = 1:nDofsPerNode
	//	            rw = nDofsPerNode*(connect(a,lmn)-1)+i;
	//	            cl = nDofsPerNode*(connect(b,lmn)-1)+k;
	//	            Stif(rw,cl) = Stif(rw,cl) + kel(nDofsPerNode*(a-1)+i,nDofsPerNode*(b-1)+k);
	//	          end
	//	        end
	//	      end
	//	    end

	//	cout << " " << endl;
	//	cout << " INSIDE ASSEMBLY" << endl;
	//	matrix.print();
	//	cout << " " << endl;

}

//-----------------------------------------------------------------
void FE_Engine::compute_rhs_vector(DENS_VEC &bc, DENS_VEC &value)
{
	int rw, cl;
	int nDofsPerNode = 2;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	int nfix = value.size();

	DENS_VEC a(nfix);
	DENS_VEC b(nfix);

	int infx = 0;
	for(int i = 0; i < (nfix/2); ++i) {
		// node number
		a(infx++) = bc(i);
		a(infx++) = bc(i);
		// dofs (0,1 or 3)
		for(int j = 0; j < nDofsPerNode; ++j) {
			b(j + infx - 2) = j;
		}
	}
	//F.reset(nTotalDofs);

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * a(n) + b(n);
		rhs( rw ) = value(n);
	}
}

void FE_Engine::compute_dirichlet_bc(DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F)
{
	int rw, cl;
	int nDofsPerNode = 2;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	int nfix = value.size();

	DENS_VEC a(nfix);
	DENS_VEC b(nfix);

	int infx = 0;
	for(int i = 0; i < (nfix/2); ++i) {
		// node number
		a(infx++) = bc(i);
		a(infx++) = bc(i);
		// dofs (0,1 or 3)
		for(int j = 0; j < nDofsPerNode; ++j) {
			b(j + infx - 2) = j;
		}
	}
	//F.reset(nTotalDofs);

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * a(n) + b(n);
		for(int cl = 0; cl < nTotalDofs; ++cl) {
			F[cl] =  ( F[cl] - GlbStiffMatrix(rw,cl) * value(n) );
		}
	}

	//cout << " COEF " << endl;
	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * a(n) + b(n);
		//		K.set( rw, rw, 1.0);

		for(int cl = 0; cl < nTotalDofs; ++cl) {
			//cout << " rw " << " " << rw << " " << " cl " << cl << endl;
			//			if (rw != cl) {
			//				K.set( rw, cl, 0.0 );
			//				K.set( cl, rw, 0.0 );
			GlbStiffMatrix( rw, cl ) = 0.0;
			GlbStiffMatrix( cl, rw ) = 0.0;
			//			}
		}
		GlbStiffMatrix( rw, rw ) = 1.0;
		//		M12.add( rw, rw, 0.0 );
		//		M12(rw)    = 1.;
		F[rw] = value(n);
	}
}


void FE_Engine::compute_dirichlet_bc_y(DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F)
{
	int rw, cl;
	int nDofsPerNode = 2;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	int nfix = value.size();

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * bc(n) + 1;
		for(int cl = 1; cl < nTotalDofs; ++cl) {
			F[cl] =  ( F[cl] - GlbStiffMatrix(rw,cl) * value(n) );
		}
	}

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * bc(n) + 1;
		for(int cl = 1; cl < nTotalDofs; ++cl) {
			GlbStiffMatrix( rw, cl ) = 0.0;
			GlbStiffMatrix( cl, rw ) = 0.0;
		}
		GlbStiffMatrix( rw, rw ) = 1.0;

		F[rw] = value(n);
	}
}


void FE_Engine::compute_dirichlet_bc_x(DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F)
{
	int rw, cl;
	int nDofsPerNode = 2;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	int nfix = value.size();

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * bc(n) + 0;
		for(int cl = 0; cl < nTotalDofs; ++cl) {
			F[cl] =  ( F[cl] - GlbStiffMatrix(rw,cl) * value(n) );
		}
	}

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * bc(n) + 0;
		for(int cl = 0; cl < nTotalDofs; ++cl) {
			GlbStiffMatrix( rw, cl ) = 0.0;
			GlbStiffMatrix( cl, rw ) = 0.0;
		}
		GlbStiffMatrix( rw, rw ) = 1.0;

		F[rw] = value(n);
	}
}


void FE_Engine::compute_dirichlet_presc(DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F, DENS_MAT &K)
{
	int rw, cl;
	int nDofsPerNode = 2;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	int nfix = value.size();

	DENS_VEC a(nfix);
	DENS_VEC b(nfix);

	int infx = 0;
	for(int i = 0; i < (nfix/2); ++i) {
		// node number
		a(infx++) = bc(i);
		a(infx++) = bc(i);
		for(int j = 0; j < nDofsPerNode; ++j) {
			b(j + infx - 2) = j;
		}
	}

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * a(n) + b(n);
		for(int cl = 0; cl < nTotalDofs; ++cl) {
			F[cl] =  ( F[cl] - K(rw,cl) * value(n) );
		}
	}

	for(int n = 0; n < nfix; ++n) {
		rw = nDofsPerNode * a(n) + b(n);

		for(int cl = 0; cl < nTotalDofs; ++cl) {
			K( rw, cl ) = 0.0;
			K( cl, rw ) = 0.0;
		}
		K( rw, rw ) = 1.0;
		F[rw] = value(n);
	}
}

void FE_Engine::compute_robin(DENS_MAT Neighbor, DENS_VEC &bc, DENS_VEC &value, DENS_VEC &F, DENS_MAT &K, double rcoeff)
{
  int rw, cl;
  int nDofsPerNode = 2;
  int nTotalDofs   = GlbStiffMatrix.nRows();

  int nfix = value.size();
  int idrch = 0;
  int FENode, NLt, NRt;

  for(int i = 0; i < bc.size(); ++i) {
    FENode = bc(i);
    NLt = Neighbor(1,FENode);
    NRt = Neighbor(2,FENode);


    K( FENode*2, FENode*2 )  += 2./3.*rcoeff; // diagonal
    K( FENode*2, NLt*2 ) += 1./6.*rcoeff;
    K( FENode*2, NRt*2 ) += 1./6.*rcoeff;

    K( FENode*2+1, FENode*2+1 )  += 2./3.*rcoeff; // diagonal
    K( FENode*2+1, NLt*2+1 ) += 1./6.*rcoeff;
    K( FENode*2+1, NRt*2+1 ) += 1./6.*rcoeff;
  }
}



void FE_Engine::FEsetup(double mu,
		double nu,
		int ndof,
		int ncoord,
		double load)
{
	time_t start,end;

	// BC application
	double dxv = -1.0;
	DENS_VEC dirichletnodes;
	DENS_VEC dirichlet_values;
	feMesh->get_Nodes_x(dxv, dirichletnodes);

	// BC application
	dxv = -1.0;
	DENS_VEC dirichletnodes1;
	DENS_VEC dirichlet_values1;
	feMesh->get_Nodes_y(dxv, dirichletnodes1);

	// Loading application
	dxv = 1.0;
	DENS_VEC rhsnodes;
	DENS_VEC rhs_values;
	feMesh->get_Nodes_x(dxv, rhsnodes);

	// -------------------------------------------------------------
	// Material obj
	// -------------------------------------------------------------
	Material Material(mu, nu, ndof, ncoord);

	// -------------------------------------------------------------
	// Quadrilateral obj
	// -------------------------------------------------------------
	FE_ElementQuad feElement;


	// -------------------------------------------------------------
	// assembly global stiffness matrix
	// -------------------------------------------------------------
	time (&start);
	fprintf ( stdout, "\n\n Computing gobal stiffness matrix... " );
	compute_stiffness_matrix(Material, feElement);
	time (&end);
	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );

	// -------------------------------------------------------------
	// compute compute dirichlet bc */
	// -------------------------------------------------------------
	dirichlet_values.reset(dirichletnodes.size());
	int idrch = 0;
	for(int i = 0; i < dirichletnodes.size(); ++i) {
		dirichlet_values(idrch++) = 0.0;
	}
	DENS_VEC Fdrchlt;
	int nTotalDofs   = GlbStiffMatrix.nRows();
	Fdrchlt.reset(nTotalDofs);
	compute_dirichlet_bc_x(dirichletnodes, dirichlet_values, Fdrchlt);

	dirichlet_values1.reset(dirichletnodes.size());
	idrch = 0;
	for(int i = 0; i < dirichletnodes.size(); ++i) {
		dirichlet_values1(idrch++) = 0.0;
	}
	DENS_VEC Fdrchlt1;
	Fdrchlt1.reset(nTotalDofs);
	compute_dirichlet_bc_y(dirichletnodes1, dirichlet_values1, Fdrchlt1);

	// -------------------------------------------------------------
	// compute right hand side */
	// -------------------------------------------------------------
	rhs_values.reset(2 * rhsnodes.size());
	int irhs = 0;
	double dfsload = load;// /rhsnodes.size();
        fprintf(stdout, "rhsnodes.size()=%d\n",rhsnodes.size());
	for(int i = 0; i < rhsnodes.size(); ++i) {
		rhs_values(irhs++) = dfsload;
		rhs_values(irhs++) = 0.0;
	}
	rhs.reset(nTotalDofs);
	compute_rhs_vector(rhsnodes, rhs_values);

	int Corner_node = feMesh->get_Node_xy(1.0, 1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(1.0, -1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

}


void FE_Engine::FEsetup1(double mu,
		double nu,
		int ndof,
		int ncoord,
		double load)
{
	time_t start,end;

	// Loading application
	double dxv = 1.0;
	DENS_VEC rhsnodes;
	DENS_VEC rhs_values;
	feMesh->get_Nodes_x(dxv, rhsnodes);

	// Loading application
	dxv = -1.0;
	DENS_VEC rhsnodes1;
	DENS_VEC rhs_values1;
	feMesh->get_Nodes_x(dxv, rhsnodes1);

	// -------------------------------------------------------------
	// Material obj
	// -------------------------------------------------------------
	Material Material(mu, nu, ndof, ncoord);

	// -------------------------------------------------------------
	// Quadrilateral obj
	// -------------------------------------------------------------
	FE_ElementQuad feElement;


	// -------------------------------------------------------------
	// assembly global stiffness matrix
	// -------------------------------------------------------------
	time (&start);
	fprintf ( stdout, "\n\n Computing gobal stiffness matrix... " );
	compute_stiffness_matrix(Material, feElement);
	time (&end);
	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );

	int nTotalDofs   = GlbStiffMatrix.nRows();

	// BC application
	DENS_VEC dirichletnodes;
	DENS_VEC dirichlet_values;
	int dirichletN = feMesh->get_Node_xy(0.0, 1.0);

	// here I'm fixing only two nodex
	dirichletnodes.reset(2);
	dirichletnodes(0) = dirichletN;
	dirichletN = feMesh->get_Node_xy(0.0, -1.0);
	dirichletnodes(1) = dirichletN;
	dirichlet_values.reset(2);
	dirichlet_values(0) = 0.0;
	dirichlet_values(1) = 0.0;

	DENS_VEC Fdrchlt;
	Fdrchlt.reset(nTotalDofs);
	compute_dirichlet_bc_x(dirichletnodes, dirichlet_values, Fdrchlt);

	dirichletN = feMesh->get_Node_xy(1.0, 0.0);
	dirichletnodes(0) = dirichletN;
	dirichletN = feMesh->get_Node_xy(-1.0, 0.0);
	dirichletnodes(1) = dirichletN;
	dirichlet_values(0) = 0.0;
	dirichlet_values(1) = 0.0;
	compute_dirichlet_bc_y(dirichletnodes, dirichlet_values, Fdrchlt);

	// -------------------------------------------------------------
	// compute right hand side */
	// -------------------------------------------------------------
	rhs.reset(nTotalDofs); // global load vector

	rhs_values.reset(2 * rhsnodes.size());
	int irhs = 0;
	double dfsload = load/rhsnodes.size();
	for(int i = 0; i < rhsnodes.size(); ++i) {
		rhs_values(irhs++) = dfsload;
		rhs_values(irhs++) = 0.0;
	}
	compute_rhs_vector(rhsnodes, rhs_values);

	rhs_values1.reset(2 * rhsnodes1.size());
	irhs = 0;
	for(int i = 0; i < rhsnodes1.size(); ++i) {
		rhs_values1(irhs++) = -dfsload;
		rhs_values1(irhs++) = 0.0;
	}
	compute_rhs_vector(rhsnodes1, rhs_values1);

	int Corner_node = feMesh->get_Node_xy(1.0, 1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(1.0, -1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(-1.0, -1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(-1.0, 1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	//	ofstream Sfile;
	//	Sfile.open ("StiffnessMat.txt");


	//	for(int i = 0; i < nTotalDofs; ++i) {
	//		for(int j = 0; j < nTotalDofs; ++j) {
	//			Sfile << " " << GlbStiffMatrix( i, j );
	//		}
	//		Sfile << " " << endl;
	//	}

	//	Sfile.close();
}



void FE_Engine::FEsetup2(double mu,
		double nu,
		int ndof,
		int ncoord,
		double load)
{
	time_t start,end;

	// Loading application (x=1,y)
	double dxv = 1.0;
	DENS_VEC rhsnodes;
	DENS_VEC rhs_values;
	feMesh->get_Nodes_x(dxv, rhsnodes);

	// Loading application (x=-1,y)
	dxv = -1.0;
	DENS_VEC rhsnodes1;
	DENS_VEC rhs_values1;
	feMesh->get_Nodes_x(dxv, rhsnodes1);
    
 	// Loading application (x,y=1)
    dxv = 1.0;
	DENS_VEC rhsnodes2;
	DENS_VEC rhs_values2;
	feMesh->get_Nodes_y(dxv, rhsnodes2);

	// Loading application (x,y=-1)
	dxv = -1.0;
	DENS_VEC rhsnodes3;
	DENS_VEC rhs_values3;
	feMesh->get_Nodes_y(dxv, rhsnodes3);
    

	// -------------------------------------------------------------
	// Material obj
	// -------------------------------------------------------------
	Material Material(mu, nu, ndof, ncoord);

	// -------------------------------------------------------------
	// Quadrilateral obj
	// -------------------------------------------------------------
	FE_ElementQuad feElement;


	// -------------------------------------------------------------
	// assembly global stiffness matrix
	// -------------------------------------------------------------
	time (&start);
	fprintf ( stdout, "\n\n Computing gobal stiffness matrix... " );
	compute_stiffness_matrix(Material, feElement);
	time (&end);
	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );

	int nTotalDofs   = GlbStiffMatrix.nRows();

	// BC application
	DENS_VEC dirichletnodes(1);
	DENS_VEC dirichlet_values(1);
	int dirichletN = feMesh->get_Node_xy(0.0, 0.0);

	// here I'm fixing only 1 nodex
//	dirichletnodes.reset(1);
	dirichletnodes(0) = dirichletN;
//	dirichletN = feMesh->get_Node_xy(0.0, -1.0);
//	dirichletnodes(1) = dirichletN;
//	dirichlet_values.reset(1);
	dirichlet_values(0) = 0.0;
//	dirichlet_values(1) = 0.0;

	DENS_VEC Fdrchlt;
	Fdrchlt.reset(nTotalDofs);
	compute_dirichlet_bc_x(dirichletnodes, dirichlet_values, Fdrchlt);

//	dirichletN = feMesh->get_Node_xy(0.0, 0.0);
//	dirichletnodes(0) = dirichletN;
//	dirichletN = feMesh->get_Node_xy(-1.0, 0.0);
//	dirichletnodes(1) = dirichletN;
//	dirichlet_values(0) = 0.0;
//	dirichlet_values(1) = 0.0;
	compute_dirichlet_bc_y(dirichletnodes, dirichlet_values, Fdrchlt);

	// BC application
    // Yue, if you look at the FEsetup1 I applied some dirichlet restriction 
    // in order to keep some simetry and to make sure that I'g get convergence
    // I don't see any simetry to the shear problem other then the central 
    // point of the geometry been fixed. The problem is that there is no FE
    // node there. I'm not sure if the FE solver will converge. I didn't
    // have time to test it.

	// -------------------------------------------------------------
	// compute right hand side */
	// -------------------------------------------------------------
	rhs.reset(nTotalDofs); // global load vector

	rhs_values.reset(2 * rhsnodes.size());
	int irhs = 0;
	double dfsload = load/rhsnodes.size();
	for(int i = 0; i < rhsnodes.size(); ++i) {
		rhs_values(irhs++) = 0.0;
		rhs_values(irhs++) = -dfsload;
	}
	compute_rhs_vector(rhsnodes, rhs_values);

	rhs_values1.reset(2 * rhsnodes1.size());
	irhs = 0;
	for(int i = 0; i < rhsnodes1.size(); ++i) {
		rhs_values1(irhs++) = 0.0;
		rhs_values1(irhs++) = dfsload;
	}
	compute_rhs_vector(rhsnodes1, rhs_values1);
    
    rhs_values2.reset(2 * rhsnodes2.size());
	irhs = 0;
	for(int i = 0; i < rhsnodes2.size(); ++i) {
		rhs_values2(irhs++) = -dfsload;
		rhs_values2(irhs++) = 0.0;
	}
	compute_rhs_vector(rhsnodes2, rhs_values2);
    
    rhs_values3.reset(2 * rhsnodes3.size());
	irhs = 0;
	for(int i = 0; i < rhsnodes3.size(); ++i) {
		rhs_values3(irhs++) = dfsload;
		rhs_values3(irhs++) = 0.0;
	}
	compute_rhs_vector(rhsnodes3, rhs_values3);

	int Corner_node = feMesh->get_Node_xy(1.0, 1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(1.0, -1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(-1.0, -1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;

	Corner_node = feMesh->get_Node_xy(-1.0, 1.0);
	rhs(2*Corner_node+0) = rhs(2*Corner_node+0)/2.0;
	rhs(2*Corner_node+1) = rhs(2*Corner_node+1)/2.0;
}

// Before Robin
//void FE_Engine::FESolve(int tsp, DENS_VEC drchlt_prsc_nds,
//DENS_VEC drchlt_prsc_vls,
//DENS_VEC &updatedCoords,
//DENS_VEC &lambda, int flag)

void FE_Engine::FESolve(int tsp, DENS_VEC drchlt_prsc_nds,
		DENS_VEC drchlt_prsc_vls,
		DENS_VEC &updatedCoords,
		DENS_VEC &lambda, int flag, double rcoeff, DENS_MAT Neighbor)
{
	time_t start,end;
	int nTotalDofs   = GlbStiffMatrix.nRows();

	// Stiffness Matrix
	DENS_MAT GlbK;

	// assemble consistent mass (nnodes X nnodes)
	GlbK.reset(nTotalDofs, nTotalDofs);

	for (int i = 0; i < nTotalDofs; ++i) {
		for (int j = 0; j < nTotalDofs; ++j) {
			GlbK(i,j) = GlbStiffMatrix(i,j);
		}
	}

	DENS_VEC Fdrchlt_prsc;
	Fdrchlt_prsc.reset(nTotalDofs);
	if (flag == 1) {
		// -------------------------------------------------------------
		// compute compute dirichlet prescribed bc */
		// -------------------------------------------------------------
		compute_dirichlet_presc(drchlt_prsc_nds, drchlt_prsc_vls, Fdrchlt_prsc, GlbK);
	}
	else if (flag == 2){
		// -------------------------------------------------------------
		// compute compute neumann bc */
		// -------------------------------------------------------------
		for(int i = 0; i < drchlt_prsc_nds.size(); ++i) {
			int NodeID = drchlt_prsc_nds(i);
			Fdrchlt_prsc(2*NodeID + 0) = drchlt_prsc_vls(2*i + 0);
			Fdrchlt_prsc(2*NodeID + 1) = drchlt_prsc_vls(2*i + 1);
		}
	}	
	else if (flag == 3){
		// -------------------------------------------------------------
		// compute compute Robin bc */
		// -------------------------------------------------------------
		compute_robin(Neighbor, drchlt_prsc_nds, drchlt_prsc_vls, Fdrchlt_prsc, GlbK, rcoeff);
		for(int i = 0; i < drchlt_prsc_nds.size(); ++i) {
			int NodeID = drchlt_prsc_nds(i);
			Fdrchlt_prsc(2*NodeID + 0) = drchlt_prsc_vls(2*i + 0);
			Fdrchlt_prsc(2*NodeID + 1) = drchlt_prsc_vls(2*i + 1);
		}
	}

	DENS_VEC rhs_total;
	rhs_total.reset(nTotalDofs);
	for (int i = 0; i < rhs.size(); ++i) {
		rhs_total(i) = rhs(i) + Fdrchlt_prsc(i);
	}

	cout << " " << endl;
	cout << " Finite Element at time step # ....:  " << tsp << endl;

	ofstream Sfile;
	Sfile.open ("RhsMat.txt");


	for(int i = 0; i < nTotalDofs; ++i) {
		Sfile << rhs_total( i ) << endl;
	}

	Sfile.close();

	// -------------------------------------------------------------
	// solve linear system */
	// -------------------------------------------------------------
	time (&start);
	fprintf ( stdout, "\n\n Solving linear system... " );
	fprintf ( stdout, "\n # of equation...: %d ", GlbStiffMatrix.nRows());

	int maxIterations = 200;
	double tolerance  = 1e-4;
	LambdaMatrixSolverCg CGSolver(maxIterations, tolerance);

	/** execute the solver */
	//	DENS_VEC lambda;
	lambda.reset( GlbStiffMatrix.nCols() );

	//	CGSolver.execute(Fdrchlt_prsc, lambda, GlbStiffMatrix);
	CGSolver.execute(rhs_total, lambda, GlbK);

	time (&end);
	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );

	// Returns New Coordinates
	DENS_MAT Coords;
	Coords =  feMesh->nodal_coordinates();

	updatedCoords.reset(GlbStiffMatrix.nRows());

	int indx = 0;
	int indw = 0;
	for (int i = 0; i < Coords.nCols(); ++i) {
		updatedCoords(indw++) = Coords(0, i) + lambda(indx++);
		updatedCoords(indw++) = Coords(1, i) + lambda(indx++);
	}
	post_processing(tsp, lambda);
}

void FE_Engine::ComputeStress(double PoissonRatio, double YoungModulus, DENS_MAT coords_xy0, DENS_MAT coords_xy, DENS_MAT Connectivity, DENS_MAT &Stress)
{
	double One_Ni, Ni, ElemArea, Ex, Ey, Yxy, Cte, Da;
	One_Ni = 1.0 - PoissonRatio;
	Ni     = PoissonRatio;
	const int nSD = coords_xy0.nRows();
	const int nne = Connectivity.nRows();
	// Stores element coordinates
	DENS_MAT xCoords0;
	xCoords0.reset(nSD, nne);
	DENS_MAT xCoords;
	xCoords.reset(nSD, nne);	
	// Stores Jacobian matrix
	DENS_MAT Jac0;
	Jac0.reset(nSD, nSD);
	DENS_MAT Jac;
	Jac.reset(nSD, nSD);
	DENS_MAT F;
	F.reset(nSD, nSD);
	DENS_MAT E;
	E.reset(nSD, nSD);
	DENS_VEC ElemStress(nne);
	DENS_VEC PatchMeasure(coords_xy0.nCols());
	Stress.reset(nne, coords_xy0.nCols());

	int inode,ielem,isd, nn, flag;
	double raio2;

	// loop over all elements
	for(ielem = 0; ielem < Connectivity.nCols(); ++ielem) {

		flag = 0;
		for (inode = 0; inode < nne; inode++) {
			// get element connectivity
			const int id = Connectivity(inode, ielem);

			// get element coordinates (Initial)
			xCoords0(0,inode) = coords_xy0(0, id);
			xCoords0(1,inode) = coords_xy0(1, id);

			raio2 = xCoords0(0,inode)*xCoords0(0,inode) + xCoords0(1,inode)*xCoords0(1,inode);
			if (raio2<(0.1*0.1)){
				flag = 0; //it controls if there is a hole or not (if flag = 1 yes, hole)
			}
			// get element coordinates (Updated)
			xCoords(0,inode)  = coords_xy(0, id);
			xCoords(1,inode)  = coords_xy(1, id);
		}
		if (flag == 0){

			// Retuns the jacobian matrix oj initial configuration
			Jac0(0,0) = xCoords0(0,0) - xCoords0(0,2);
			Jac0(0,1) = xCoords0(1,0) - xCoords0(1,2);
			Jac0(1,0) = xCoords0(0,1) - xCoords0(0,2);
			Jac0(1,1) = xCoords0(1,1) - xCoords0(1,2);

			// Returns the Jacobian matrix compute on the current configuration
			Jac(0,0) = xCoords(0,0) - xCoords(0,2);
			Jac(0,1) = xCoords(1,0) - xCoords(1,2);
			Jac(1,0) = xCoords(0,1) - xCoords(0,2);
			Jac(1,1) = xCoords(1,1) - xCoords(1,2);

			ElemArea  = 0.5 *(xCoords0(0,0)*xCoords0(1,1) + xCoords0(0,2)*xCoords0(1,0) + xCoords0(0,1)*xCoords0(1,2));
			ElemArea -= 0.5 *(xCoords0(0,2)*xCoords0(1,1) + xCoords0(0,0)*xCoords0(1,2) + xCoords0(0,1)*xCoords0(1,0));

			//F = (JacobianMatrix0 \ JacobianMatrix1)';
			Da  = 1.0 / (  Jac0(0,0)*Jac0(1,1)-Jac0(0,1)*Jac0(1,0));
			F(0,0) = Da * ( Jac(0,0)*Jac0(1,1)-Jac0(0,1)* Jac(1,0));
			F(0,1) = Da * (Jac0(0,0)* Jac(1,0)- Jac(0,0)*Jac0(1,0));
			F(1,0) = Da * ( Jac(0,1)*Jac0(1,1)-Jac0(0,1)* Jac(1,1));
			F(1,1) = Da * (Jac0(0,0)* Jac(1,1)- Jac(0,1)*Jac0(1,0));

			// Get Green-Lagrange strain tensor
			// E = 0.5*(F + F') - eye(size(F));
			E(0,0) = F(0,0)-1.0;
			E(0,1) = 0.5*(F(0,1)+F(1,0));
			E(1,0) = E(0,1);
			E(1,1) = F(1,1)-1.0;

			// 'PLANE_STRAIN'
			// Get Green-Lagrange strain tensor components
			Ex  = E(0,0);
			Yxy = E(0,1) + E(1,0);
			Ey  = E(1,1);

			Cte = YoungModulus / ((1.0 + Ni) * (One_Ni - Ni));

			ElemStress(0) = Cte * (One_Ni * Ex + Ni * Ey);     // Sx
			ElemStress(1) = Cte * (Ni * Ex + One_Ni * Ey);     // Sy
			ElemStress(2) = Cte * ((One_Ni - Ni) / 2.0) * Yxy; // Sxy

			for (inode = 0; inode < nne; inode++) {
				// get element connectivity
				const int id = Connectivity(inode, ielem);						
				PatchMeasure( id ) = PatchMeasure( id )+ElemArea;			
				Stress(0, id) = Stress(0, id)+ElemStress(0)*ElemArea;
				Stress(1, id) = Stress(1, id)+ElemStress(1)*ElemArea;
				Stress(2, id) = Stress(2, id)+ElemStress(2)*ElemArea;
			}
		}
	}

	// Nodal stress
	for (nn = 0; nn < coords_xy0.nCols(); nn++) {
		for (inode = 0; inode < nne; inode++) {
			Stress(inode,nn)  = Stress(inode,nn) / PatchMeasure(nn);
		}
	}

	//cout << "PRINT STRESS with ELMTarea:" << ElemArea << " PatchMeasure:" << PatchMeasure(0) << " Jac0(0,0)=" << Jac0(0,0) << " Jac(0,0)=" << Jac(0,0) << endl;
	//    cout << "YM=" << YoungModulus << endl;
	//for(int i = 0; i < Stress.nCols()/10. ; ++i)
	//	cout << Stress(0, i) << " " << Stress(1, i) << " " << Stress(2, i) << " " << endl;

}

//* Here it's used the Delaunay to create a mesh
void FE_Engine::GenerateMesh(DENS_MAT coords_xy, DENS_MAT &Connectivity)
{
	Delaunay_Triangulation(coords_xy, Connectivity);
}


void FE_Engine::post_processing(int tsp, DENS_VEC &lambda)
{
	// -------------------------------------------------------------
	// Post-Processing - gather information for lammps*/
	// -------------------------------------------------------------
	DENS_MAT Coords, Coords2;
	Coords =  feMesh->nodal_coordinates();
	int indx = 0;
	int indw = 0;

	// -------------------------------------------------------------
	// Post-Processing */
	// -------------------------------------------------------------
	Coords2 = Coords;
	indx = 0;
	indw = 1;
	for (int i = 0; i < Coords.nCols(); ++i) {
		Coords2(0, i) = Coords(0, i) + lambda(indx++);
		Coords2(1, i) = Coords(1, i) + lambda(indw++);
		indx++;
		indw++;
	}

	DENS_MAT xyz( (Coords.nRows()+1),  Coords.nCols());
	for (int i = 0; i < Coords.nRows(); ++i) {
		for (int j = 0; j < Coords.nCols(); ++j) {
			xyz(i, j) = Coords2(i, j);
		}
	}

	Array2D<int> Connect;
	Connect = feMesh->connectivity();

	DENS_VEC p( Coords.nCols());
	DENS_MAT uvw( 3, Coords.nCols());
	int ind=0;
	for (int i = 0; i < Coords.nCols(); ++i) {
		p(i)  = lambda(ind++);
		ind++;
	}
	ind=0;
	for (int i = 0; i < Coords.nCols(); ++i) {
		uvw(0, i) = lambda(ind++);
		uvw(1, i) = lambda(ind++);
		uvw(2, i) = 0.0;
	}

	twod_to_vtk (tsp, xyz, Connect, p, uvw);
	//	twod_to_vtk (DENS_MAT &xyz, Array2D<int> &e_conn, DENS_VEC &p, DENS_MAT &uvw);
}


void FE_Engine::twod_to_vtk (int tsp, DENS_MAT &xyz, Array2D<int> &e_conn, DENS_VEC &p, DENS_MAT &uvw)
//void FE_Engine::twod_to_vtk ( xy, e_conn, p, uvw, output_filename, title )
//*****************************************************************************80
//
//// TWOD_TO_VTK writes out a TWOD dataset to a legacy VTK file.
//
//  Discussion:
//
//    The VTK file can be read and displayed by the Paraview program.
//
//    Thanks to Mike Sussman for suggesting that real data should be
//    written with the "double" attribute rather than the "float",
//    JVB, 20 December 2010.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real XY(NODE_NUM,2), the node coordinates.
//
//    Input, integer E_CONN(ELEMENT_NUM,ELEMENT_ORDER), the nodes that
//    form each element.  Node indices are 1-based.
//
//    Input, real UVP(NODE_NUM,3), U and V velocity components and
//    pressure at each node.
//
//    Input, string OUTPUT_FILENAME, the name of the output file.
//    By convention, this file should have the extension ".vtk".
//
//    Input, string TITLE, a title for the data.
//

//
//  Determine the sizes of things.
//
{
	FILE* output_unit;

	int nNodes = feMesh->get_nNodes();
	int dim_num = 2;

	int nElements = feMesh->get_nElements();
	int element_order =  4;
	//
	//  Open the output file.
	//
	char str[80];
	std::string s = number_to_string(tsp);
	//	std::string s = number_to_string(0);

	char *stsp =new char[s.size()+1];
	stsp[s.size()]=0;
	memcpy(stsp,s.c_str(),s.size());

	strcpy (str,"ns2d_fem_");
	strcat (str, stsp);
	strcat (str, ".vtk");


	int Choose = 2;
	// Prints solution.
	if (Choose == 1) {
		output_unit = stdout;
	} else if ((output_unit = fopen(str, "w")) == NULL) {
		cout << "Error in opening %s file, ns2d_fem.vtk" << endl;
	}

	if ( element_order == 6 ) {
		fprintf ( stdout, "\n" );
		fprintf ( stdout, "TWO_TO_VTK:\n" );
		fprintf ( stdout, "  The input data uses quadratic elements.\n" );
		fprintf ( stdout, "  The output data will use linear elements.\n" );
	}

	//
	//  Write the data.
	//
	//  vtk_puv_write (output_unit, title, nNodes, nElements, ...
	//    element_order2, xyz, e_conn, p, uvw );

	vtk_puv_write (output_unit, nNodes, nElements, element_order, xyz, e_conn, p, uvw );

	fclose ( output_unit );

	fprintf ( stdout, "\n The data was written to *.vtk" );
	cout << " " << endl;

}

//void FE_Engine::vtk_puv_write ( output_unit, title, nNodes, nElements, element_order, xyz, element_node, p, uvw )
void FE_Engine::vtk_puv_write (FILE* output_unit, int nNodes, int nElements, int element_order, DENS_MAT &xyz, Array2D<int> &e_conn, DENS_VEC &p, DENS_MAT &uvw ) const

//*****************************************************************************80
//
//// VTK_PUV_WRITE writes pressure and velocity data to a VTK file.
//
//  Discussion:
//
//    The data is assumed to have been computed by a finite element program
//    for a 2D geometry which has been meshed using triangular elements
//    of 3 or 6 nodes.
//
//    The solution data includes the pressure and velocity vector at each node.
//
//    Note that the VTK format used here is known as the "legacy" or "old style"
//    format.  It has been superseded by a family of XML based formats.  The
//    appropriate replacement for the VTK format used here is known as "VTU",
//    which is the Visual Toolkit format for unstructured grid data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer OUTPUT_UNIT, the output unit.
//
//    Input, string TITLE, a title for the data.
//
//    Input, integer NODE_NUM, the number of nodes.
//
//    Input, integer ELEMENT_NUM, the number of elements.
//
//    Input, integer ELEMENT_ORDER, the order of the elements.
//
//    Input, real XYZ(3,NODE_NUM), the node coordinates.
//
//    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
//    nodes that make up each element.  Node indices are zero-based.
//
//    Input, real P(1,NODE_NUM), the pressure at each node.
//
//    Input, real UVW(3,NODE_NUM), the velocity at each node.
//

{
	fprintf ( output_unit, "# vtk DataFile Version 2.0\n" );
	fprintf ( output_unit, "UVP for 10x1 bending problem \n");
	fprintf ( output_unit, "ASCII\n" );
	fprintf ( output_unit, "\n" );
	fprintf ( output_unit, "DATASET UNSTRUCTURED_GRID\n" );
	fprintf ( output_unit, "POINTS %d double\n", nNodes );

	for(int node = 0; node < nNodes; ++node) {
		fprintf ( output_unit, "  %f  %f  %f\n", xyz(0,node), xyz(1,node), xyz(2,node) );
	}
	//
	//  Note that CELL_SIZE uses ELEMENT_ORDER+1 because the order of each element
	//  is included as a data item.
	//
	int cell_size = nElements * ( element_order + 1 );

	fprintf ( output_unit, "\n" );
	fprintf ( output_unit, "CELLS  %d  %d\n", nElements, cell_size );
	for(int element = 0; element < nElements; ++element) {
		fprintf ( output_unit, "  %d", element_order );
		for(int order = 0; order < element_order; ++order) {
			fprintf ( output_unit, "  %d", e_conn(order, element) );
		}
		fprintf ( output_unit, "\n" );
	}

	//
	// VTK_ QUAD (= 9)
	//
	fprintf ( output_unit, "\n" );
	fprintf ( output_unit, "CELL_TYPES %d\n", nElements );

	for(int element = 0; element < nElements; ++element) {
		fprintf ( output_unit, "9\n" );
	}


	fprintf ( output_unit, "\n" );
	fprintf ( output_unit, "POINT_DATA %d\n", nNodes );
	//	fprintf ( output_unit, "SCALARS pressure double\n" );
	//	fprintf ( output_unit, "LOOKUP_TABLE default\n" );
	//	for(int node = 0; node < nNodes; ++node) {
	//		fprintf ( output_unit, "  %f\n", p(node) );
	//	}
	fprintf ( output_unit, "VECTORS displacement double\n" );
	for(int node = 0; node < nNodes; ++node)  {
		fprintf ( output_unit, "  %f  %f  %f\n", uvw(0, node), uvw(1, node), uvw(2, node) );
	}
}



template <typename T>
std::string FE_Engine::number_to_string(T number)
{
	return dynamic_cast<std::stringstream *> (&(std::stringstream() << number))->str();
}

void FE_Engine::compute_neumann_bc() const
{

}

void FE_Engine::driver(	double mu,
		double nu,
		int ndof,
		int ncoord,
		double load,
		DENS_VEC drchlt_prsc_nds,
		DENS_VEC drchlt_prsc_vls,
		DENS_VEC &updatedCoords,
		DENS_VEC &lambda) const
{
	//	time_t start,end;
	////	// -------------------------------------------------------------
	////	// Mesh obj
	////	// -------------------------------------------------------------
	////	time (&start);
	////	fprintf ( stdout, "\n Generating mesh... " );
	////	FE_Uniform2DMesh feMesh(nx_, ny_, xmin_, xmax_, ymin_, ymax_, xscale_, yscale_);
	////	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//
	//	time (&start);
	//	double dxv = -1.0;
	//	DENS_VEC dirichletnodes;
	//	DENS_VEC dirichlet_values;
	////	fprintf ( stdout, "\n\n Finding nodes with Dirichlet BC... " );
	//	feMesh->get_Nodes_x(dxv, dirichletnodes);
	////	dirichletnodes.print();
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//
	//	dxv = -1.0;
	//	DENS_VEC dirichletnodes1;
	//	DENS_VEC dirichlet_values1;
	//	feMesh->get_Nodes_y(dxv, dirichletnodes1);
	//
	//
	//	time (&start);
	//	dxv = 1;
	//	DENS_VEC rhsnodes;
	//	DENS_VEC rhs_values;
	////	fprintf ( stdout, "\n\n Finding loded nodes... " );
	//	feMesh->get_Nodes_x(dxv, rhsnodes);
	////	rhsnodes.print();
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//
	//	// -------------------------------------------------------------
	//	// Material obj
	//	// -------------------------------------------------------------
	//	Material Material(mu, nu, ndof, ncoord);
	//
	//
	//	// -------------------------------------------------------------
	//	// Quadrilateral obj
	//	// -------------------------------------------------------------
	//	FE_ElementQuad feElement;
	//
	//
	//	// -------------------------------------------------------------
	//	// assembly global stiffness matrix
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Computing gobal stiffness matrix... " );
	////	SPAR_MAT GlbStiffMatrix;
	//	DENS_MAT GlbStiffMatrix;
	//	compute_stiffness_matrix(Material, feElement, GlbStiffMatrix);
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//	// -------------------------------------------------------------
	//	// compute compute dirichlet bc */
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Applying Dirichlet BC... " );
	//	dirichlet_values.reset(ndof * dirichletnodes.size());
	//	int idrch = 0;
	//	for(int i = 0; i < dirichletnodes.size(); ++i) {
	//		dirichlet_values(idrch++) = 0.0;
	//		dirichlet_values(idrch++) = 0.0;
	//	}
	//	DENS_VEC Fdrchlt;
	//	int nTotalDofs   = GlbStiffMatrix.nRows();
	//	Fdrchlt.reset(nTotalDofs);
	//	compute_dirichlet_bc(dirichletnodes, dirichlet_values, Fdrchlt, GlbStiffMatrix);
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//
	//	dirichlet_values1.reset(dirichletnodes.size());
	//	idrch = 0;
	//	for(int i = 0; i < dirichletnodes.size(); ++i) {
	//		dirichlet_values1(idrch++) = 0.0;
	//	}
	//	DENS_VEC Fdrchlt1;
	//	Fdrchlt1.reset(nTotalDofs);
	//	compute_dirichlet_bc_y(dirichletnodes1, dirichlet_values1, Fdrchlt1, GlbStiffMatrix);
	//
	//
	//	// -------------------------------------------------------------
	//	// compute compute dirichlet prescribed bc */
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Applying Dirichlet BC... " );
	//	DENS_VEC Fdrchlt_prsc;
	//	Fdrchlt_prsc.reset(nTotalDofs);
	//	compute_dirichlet_bc(drchlt_prsc_nds, drchlt_prsc_vls, Fdrchlt_prsc, GlbStiffMatrix);
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//
	//	// -------------------------------------------------------------
	//	// compute right hand side */
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Computing right hand side... " );
	//	rhs_values.reset(2 * rhsnodes.size());
	//	int irhs = 0;
	//	double dfsload = load/rhsnodes.size();
	//	for(int i = 0; i < rhsnodes.size(); ++i) {
	//		rhs_values(irhs++) = dfsload;
	//		rhs_values(irhs++) = 0.0;
	//	}
	//	DENS_VEC rhs;
	//	rhs.reset(nTotalDofs);
	//	compute_rhs_vector(rhsnodes, rhs_values, rhs, GlbStiffMatrix);
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//	for (int i = 0; i < rhs.size(); ++i) {
	////		rhs(i) = rhs(i) + rhs_values(i);
	//		rhs(i) = rhs(i) + Fdrchlt_prsc(i);
	////		rhs(i) = rhs(i) + Fdrchlt1(i);
	////		rhs(i) = rhs(i) + Fdrchlt(i);
	//	}
	//
	//	// -------------------------------------------------------------
	//	// solve linear system */
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Solving linear system... " );
	////	fprintf ( stdout, "\n # of equation...: %d ", GlbStiffMatrix.nRows());
	//
	//	int maxIterations = 200;
	//	double tolerance  = 1000e-10;
	//	LambdaMatrixSolverCg CGSolver(maxIterations, tolerance);
	//
	//	/** execute the solver */
	////	DENS_VEC lambda;
	//	lambda.reset( GlbStiffMatrix.nCols() );
	////	CGSolver.execute(Fdrchlt_prsc, lambda, GlbStiffMatrix);
	//	CGSolver.execute(rhs, lambda, GlbStiffMatrix);
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	//
	//	// -------------------------------------------------------------
	//	// Post-Processing - gather information for lammps*/
	//	// -------------------------------------------------------------
	//	DENS_MAT Coords, Coords2;
	//	Coords =  feMesh->nodal_coordinates();
	//	updatedCoords.reset(GlbStiffMatrix.nRows());
	//	int indx = 0;
	//	int indw = 0;
	//	for (int i = 0; i < Coords.nCols(); ++i) {
	//		updatedCoords(indw++) = Coords(0, i) + lambda(indx++);
	//		updatedCoords(indw++) = Coords(1, i) + lambda(indx++);
	////		updatedCoords(indw++) = lambda(indx++);
	////		updatedCoords(indw++) = lambda(indx++);
	//	}
	//
	//	// -------------------------------------------------------------
	//	// Post-Processing */
	//	// -------------------------------------------------------------
	//	time (&start);
	////	fprintf ( stdout, "\n\n Generating vtk output file... " );
	//	Coords2 = Coords;
	//	indx = 0;
	//	indw = 1;
	//	for (int i = 0; i < Coords.nCols(); ++i) {
	//		Coords2(0, i) = Coords(0, i) + lambda(indx++);
	//		Coords2(1, i) = Coords(1, i) + lambda(indw++);
	//		indx++;
	//		indw++;
	//	}
	//
	//	DENS_MAT xyz( (Coords.nRows()+1),  Coords.nCols());
	//	for (int i = 0; i < Coords.nRows(); ++i) {
	//		for (int j = 0; j < Coords.nCols(); ++j) {
	//			xyz(i, j) = Coords2(i, j);
	//		}
	//	}
	//
	//	Array2D<int> Connect;
	//	Connect = feMesh->connectivity();
	//
	//	DENS_VEC p( Coords.nCols());
	//	DENS_MAT uvw( 3, Coords.nCols());
	//	int ind=0;
	//	for (int i = 0; i < Coords.nCols(); ++i) {
	//		p(i)  = lambda(ind++);
	//		ind++;
	//	}
	//	ind=0;
	//	for (int i = 0; i < Coords.nCols(); ++i) {
	//		uvw(0, i) = lambda(ind++);
	//		uvw(1, i) = lambda(ind++);
	//		uvw(2, i) = 0.0;
	//	}
	//
	//	twod_to_vtk (xyz, Connect, p, uvw);
	////	twod_to_vtk (DENS_MAT &xyz, Array2D<int> &e_conn, DENS_VEC &p, DENS_MAT &uvw);
	//
	//	time (&end);
	////	fprintf ( stdout, "\n Time elapsed...:  %f (sec) ", difftime (end,start) );
	////
	//	fprintf ( stdout, "\n\n Ending FE execution... \n\n" );
}


// FROM HERE TABLE_DELAUNAY
//****************************************************************************80

void Compute_Pressure(int ncoord, int ndof, int nfacenodes, int TotalNumElem, int * GlobalNodeNumber, int * connect, double * coords, double * traction, double ** xi, double * wi, double ** GlobalTraction)

//****************************************************************************80
//
// Compute the pressure on a given a region incidence
//
{
	int npoints = 3;
	int rw;
	double dt;

	double ** dxdxi = new double*[ncoord];
	for (int i = 0; i < ncoord; ++i) dxdxi[i] = new double[ncoord];
	for (int i = 0; i < ncoord; ++i) {
		for (int j = 0; j < ncoord; ++j)
		{
			dxdxi[i][j] = 0.;
		}
	}
	double ** dxidx = new double*[ncoord];
	for (int i = 0; i < ncoord; ++i) dxidx[i] = new double[ncoord];

	double ** r = new double*[nfacenodes];
	for (int i = 0; i < nfacenodes; ++i) r[i] = new double[ndof];
	for (int i = 0; i < nfacenodes; ++i) {
		for (int j = 0; j < ndof; ++j)
		{
			r[i][j] = 0.;
		}
	}

	double ** dNdxi = new double*[nfacenodes];
	for (int i = 0; i < nfacenodes; ++i) dNdxi[i] = new double[2];


	double * N = new double[nfacenodes];

	for (int Elem = 0; Elem < TotalNumElem; ++Elem) {

		for (int intpt = 0; intpt < npoints; ++intpt) {

			shapefunctions(intpt, xi, N);
			shapefunctionderivs(dNdxi);

			//
			//     Compute the jacobian matrix && its determinant
			//
			for (int i = 0; i < ncoord; ++i) {
				for (int j = 0; j < ncoord; ++j) {
					dxdxi[i][j] = 0.;
					for (int a = 0; a <  nfacenodes; ++a) {
						dxdxi[i][j] = dxdxi[i][j] + coords[ 2*(connect[a + Elem*ndof] - 1) + i] * dNdxi[a][j];
					}
				}
			}

			double dt1 = 1./(dxdxi[0][0] * dxdxi[1][1] - dxdxi[0][1] * dxdxi[1][0]);

			double dt2 = dxdxi[1][1];
			double dt3 = dxdxi[0][0];
			double dt4 = dxdxi[0][1];
			double dt5 = dxdxi[1][0];

			dxidx[0][0] = dt2/dt1;
			dxidx[1][1] = dt3/dt1;
			dxidx[0][1] = -dt4/dt1;
			dxidx[1][0] = -dt5/dt1;

			dt = sqrt( dxidx[0][0] * dxidx[1][1] - dxidx[0][1] * dxidx[1][0]);
			cout << "Jacobiano  " << dt << "\n";

			for (int a = 0; a < nfacenodes; ++a) {
				r[a][0] = r[a][0] + N[a]*traction[0]*wi[intpt]*dt;
				r[a][1] = r[a][1] + N[a]*traction[1]*wi[intpt]*dt;
				r[a][2] = r[a][2] + N[a]*traction[2]*wi[intpt]*dt;
			}

			for (int i = 0; i < ncoord; ++i) {
				for (int j = 0; j < ncoord; ++j)
				{
					dxdxi[i][j] = 0.;
				}
			}
		}


		//
		//    Assemble the element load vector into global vector
		//
		for (int a = 0; a < nfacenodes; ++a) {
			rw = GlobalNodeNumber[ connect[a + Elem*ndof] - 1 ];
			GlobalTraction[rw][0] = GlobalTraction[rw][0] + r[a][0];
			GlobalTraction[rw][1] = GlobalTraction[rw][1] + r[a][1];
			GlobalTraction[rw][2] = GlobalTraction[rw][2] + r[a][2];
		}

		//
		// vanish (2d Arrays)
		//
		for (int i = 0; i < nfacenodes; ++i) {
			for (int j = 0; j < ndof; ++j)
			{
				r[i][j] = 0.;
			}
		}
	}

	//
	// de-allocate memory (1D Array)
	//
	delete [] N;

	//
	// de-allocate memory (2d Arrays)
	//
	for (int i = 0; i < ncoord; ++i)
		delete [] dxdxi[i];
	delete [] dxdxi;

	for (int i = 0; i < nfacenodes; ++i)
		delete [] r[i];
	delete [] r;

	for (int i = 0; i < nfacenodes; ++i)
		delete [] dNdxi[i];
	delete [] dNdxi;

}




//****************************************************************************80

void shapefunctions(int nip, double ** xi_, double * N)

//****************************************************************************80
//
// Returns the Shape Functions for a given point
//
{
	// Shape Functions
	N[0] = xi_[0][nip];
	N[1] = xi_[1][nip];
	N[2] = 1.-xi_[0][nip]-xi_[1][nip];


}



//****************************************************************************80

void shapefunctionderivs(double ** dNdxi)

//****************************************************************************80
//
// Returns the Shape Function derivatives for a given point
//
{
	dNdxi[0][0] =  1.;
	dNdxi[0][1] =  0.;
	dNdxi[1][0] =  0.;
	dNdxi[1][1] =  1.;
	dNdxi[2][0] = -1.;
	dNdxi[2][1] = -1.;
}

//****************************************************************************80

void Delaunay_Triangulation (DENS_MAT coords_xy, DENS_MAT &Connectivity)

//****************************************************************************80
//
//  Purpose:
//
//    Delaunay_Triangulation is the main program for TABLE_DELAUNAY.
//
//  Discussion:
//
//    TABLE_DELAUNAY computes the Delaunay triangulation of a TABLE dataset.
//
//    The dataset is simply a set of points in the plane.
//
//    Thus, given a set of points V1, V2, ..., VN, we apply a standard 
//    Delaunay triangulation.  The Delaunay triangulation is an organization 
//    of the data into triples, forming a triangulation of the data, with
//    the property that the circumcircle of each triangle never contains
//    another data point.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    February 20/2012 (Fabiano F. Bargos)
//
//  Author:
//
//    John Burkardt
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/table_delaunay/table_delaunay.html
//
{ 

	int node_dim = coords_xy.nRows();
	int node_num = coords_xy.nCols();

	double *node_xy = new double[node_dim * node_num];
	int offset2d = 0;
	for (int i = 0; i < node_num; ++i) {
		node_xy[offset2d + 0] = coords_xy(0,i);
		node_xy[offset2d + 1] = coords_xy(1,i);
		offset2d += 2;
	}

	int triangle_num;


	int base = 1;
	int *triangle_node;
	int *triangle_neighbor;
	int triangle_order;

	//cout << "\n";
	timestamp ( );

	//cout << "\n";
	//cout << "TABLE_DELAUNAY\n";
	//cout << "  Compute the Delaunay triangulation.\n";
	//cout << "  Write an integer TABLE dataset of the triangulation.\n";

	//  cout << "  Node dimension NODE_DIM = " << node_dim << "\n";
	//cout << "  Node number    NODE_NUM = " << node_num << "\n";

	if ( node_dim != 2 )
	{
		cout << "\n";
		cout << "TABLE_DELAUNAY - Fatal error!\n";
		cout << "  The node dimension is not 2.\n";
		exit ( 1 );
	}

	//  r8mat_transpose_print_some ( node_dim, node_num, node_xy, 1, 1, node_dim, 5, 
	//    "  Initial portion of node data:" );
	//
	//  Determine the Delaunay triangulation.
	//
	triangle_order = 3;

	triangle_node     = new int[triangle_order*3*node_num];
	triangle_neighbor = new int[triangle_order*3*node_num];

	dtris2 ( node_num, base, node_xy, &triangle_num, triangle_node, triangle_neighbor );
	//
	//  Print a portion of the triangulation.
	//
	cout << "\n";
	cout << "  Computed the triangulation.\n";
	cout << "  Number of triangles is " << triangle_num << "\n";

	//  i4mat_transpose_print_some ( triangle_order, *triangle_num, triangle_node,
	//    1, 1, 3, 5, "  Initial portion of triangulation data:" );
	//
	//  Terminate execution.
	//
	//  cout << "  TABLE_DELAUNAY:\n";
	//  cout << "  Normal end of execution.\n";
	//  timestamp ( );

	// Copy connectivity
	Connectivity.reset(triangle_order, triangle_num);

	// the -1 is because triangle_node starts in one
	offset2d = 0;
	for (int i = 0; i < triangle_num; ++i) {
		Connectivity(0, i) = triangle_node[offset2d + 0] - 1;
		Connectivity(1, i) = triangle_node[offset2d + 1] - 1;
		Connectivity(2, i) = triangle_node[offset2d + 2] - 1;
		offset2d += 3;
	}

	delete [] triangle_neighbor;
	delete [] triangle_node;
	delete [] node_xy;
}


//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
	if ( 97 <= c1 && c1 <= 122 ) 
	{
		c1 = c1 - 32;
	} 
	if ( 97 <= c2 && c2 <= 122 ) 
	{
		c2 = c2 - 32;
	}     

	return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
	int digit;

	if ( '0' <= c && c <= '9' )
	{
		digit = c - '0';
	}
	else if ( c == ' ' )
	{
		digit = 0;
	}
	else
	{
		digit = -1;
	}

	return digit;
}
//****************************************************************************80

int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2, 
		double x3, double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the 
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
	double ca;
	double cb;
	double dx10;
	double dx12;
	double dx30;
	double dx32;
	double dy10;
	double dy12;
	double dy30;
	double dy32;
	double s;
	double tol;
	double tola;
	double tolb;
	int value;

	tol = 100.0 * r8_epsilon ( );

	dx10 = x1 - x0;
	dy10 = y1 - y0;
	dx12 = x1 - x2;
	dy12 = y1 - y2;
	dx30 = x3 - x0;
	dy30 = y3 - y0;
	dx32 = x3 - x2;
	dy32 = y3 - y2;

	tola = tol * r8_max ( fabs ( dx10 ), 
			r8_max ( fabs ( dy10 ), 
					r8_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

	tolb = tol * r8_max ( fabs ( dx12 ), 
			r8_max ( fabs ( dy12 ), 
					r8_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

	ca = dx10 * dx30 + dy10 * dy30;
	cb = dx12 * dx32 + dy12 * dy32;

	if ( tola < ca && tolb < cb )
	{
		value = -1;
	}
	else if ( ca < -tola && cb < -tolb )
	{
		value = 1;
	}
	else
	{
		tola = r8_max ( tola, tolb );
		s = ( dx10 * dy30 - dx30 * dy10 ) * cb 
				+ ( dx32 * dy12 - dx12 * dy32 ) * ca;

		if ( tola < s )
		{
			value = -1;
		}
		else if ( s < -tola )
		{
			value = 1;
		}
		else
		{
			value = 0;
		}

	}

	return value;
}
//****************************************************************************80

int dtris2 ( int point_num, int base, double point_xy[], int *tri_num, 
		int tri_vert[], int tri_nabe[] )

//****************************************************************************80
//
//  Purpose:
//
//    DTRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of vertices.
//
//    Input, int BASE, the base for the indexing of TRI_VERT.
//    0, use 0-based indexing.
//    1, use 1-based indexing.
//
//    Input/output, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//    On output, the vertices have been sorted into dictionary order.
//
//    Output, int *TRI_NUM, the number of triangles in the triangulation;
//    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRI_VERT[TRI_NUM*3], the nodes that make up each triangle.
//    The elements are indices of POINT_XY.  The vertices of the triangles are
//    in counter clockwise order.
//
//    Output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRI_NABE[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int DTRIS2, is 0 for no error.
{
	double cmax;
	int e;
	int error;
	int i;
	int *indx;
	int j;
	int k;
	int l;
	int ledg;
	int lr;
	int ltri;
	int m;
	int m1;
	int m2;
	int n;
	int redg;
	int rtri;
	int *stack;
	int t;
	double tol;
	int top;

	stack = new int[point_num];

	tol = 100.0 * r8_epsilon ( );
	//
	//  Sort the vertices by increasing (x,y).
	//
	indx = r82vec_sort_heap_index_a ( point_num, base, point_xy );

	r82vec_permute ( point_num, indx, base, point_xy );
	//
	//  Make sure that the data points are "reasonably" distinct.
	//
	m1 = 1;

	for ( i = 2; i <= point_num; i++ )
	{
		m = m1;
		m1 = i;

		k = -1;

		for ( j = 0; j <= 1; j++ )
		{
			cmax = r8_max ( fabs ( point_xy[2*(m-1)+j] ), 
					fabs ( point_xy[2*(m1-1)+j] ) );

			if ( tol * ( cmax + 1.0 ) 
					< fabs ( point_xy[2*(m-1)+j] - point_xy[2*(m1-1)+j] ) ) 
			{
				k = j;
				break;
			}

		}

		if ( k == -1 )
		{
			cout << "\n";
			cout << "DTRIS2 - Fatal error!\n";
			cout << "  Fails for point number I = " << i << "\n";
			cout << "  M =  " << m  << "\n";
			cout << "  M1 = " << m1 << "\n";
			cout << "  X,Y(M)  = " << point_xy[2*(m-1)+0] << "  "
					<< point_xy[2*(m-1)+1] << "\n";
			cout << "  X,Y(M1) = " << point_xy[2*(m1-1)+0] << "  "
					<< point_xy[2*(m1-1)+1] << "\n";
			delete [] stack;
			return 224;
		}

	}
	//
	//  Starting from points M1 and M2, search for a third point M that
	//  makes a "healthy" triangle (M1,M2,M)
	//
	m1 = 1;
	m2 = 2;
	j = 3;

	for ( ; ; )
	{
		if ( point_num < j )
		{
			cout << "\n";
			cout << "DTRIS2 - Fatal error!\n";
			delete [] stack;
			return 225;
		}

		m = j;

		lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1],
				point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1], 
				point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

		if ( lr != 0 )
		{
			break;
		}

		j = j + 1;

	}
	//
	//  Set up the triangle information for (M1,M2,M), and for any other
	//  triangles you created because points were collinear with M1, M2.
	//
	*tri_num = j - 2;

	if ( lr == -1 )
	{
		tri_vert[3*0+0] = m1;
		tri_vert[3*0+1] = m2;
		tri_vert[3*0+2] = m;
		tri_nabe[3*0+2] = -3;

		for ( i = 2; i <= *tri_num; i++ )
		{
			m1 = m2;
			m2 = i+1;
			tri_vert[3*(i-1)+0] = m1;
			tri_vert[3*(i-1)+1] = m2;
			tri_vert[3*(i-1)+2] = m;
			tri_nabe[3*(i-1)+0] = -3 * i;
			tri_nabe[3*(i-1)+1] = i;
			tri_nabe[3*(i-1)+2] = i - 1;

		}

		tri_nabe[3*(*tri_num-1)+0] = -3 * (*tri_num) - 1;
		tri_nabe[3*(*tri_num-1)+1] = -5;
		ledg = 2;
		ltri = *tri_num;
	}
	else
	{
		tri_vert[3*0+0] = m2;
		tri_vert[3*0+1] = m1;
		tri_vert[3*0+2] = m;
		tri_nabe[3*0+0] = -4;

		for ( i = 2; i <= *tri_num; i++ )
		{
			m1 = m2;
			m2 = i+1;
			tri_vert[3*(i-1)+0] = m2;
			tri_vert[3*(i-1)+1] = m1;
			tri_vert[3*(i-1)+2] = m;
			tri_nabe[3*(i-2)+2] = i;
			tri_nabe[3*(i-1)+0] = -3 * i - 3;
			tri_nabe[3*(i-1)+1] = i - 1;
		}

		tri_nabe[3*(*tri_num-1)+2] = -3 * (*tri_num);
		tri_nabe[3*0+1] = -3 * (*tri_num) - 2;
		ledg = 2;
		ltri = 1;
	}
	//
	//  Insert the vertices one at a time from outside the convex hull,
	//  determine visible boundary edges, and apply diagonal edge swaps until
	//  Delaunay triangulation of vertices (so far) is obtained.
	//
	top = 0;

	for ( i = j+1; i <= point_num; i++ )
	{
		m = i;
		m1 = tri_vert[3*(ltri-1)+ledg-1];

		if ( ledg <= 2 )
		{
			m2 = tri_vert[3*(ltri-1)+ledg];
		}
		else
		{
			m2 = tri_vert[3*(ltri-1)+0];
		}

		lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1], 
				point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1], 
				point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

		if ( 0 < lr )
		{
			rtri = ltri;
			redg = ledg;
			ltri = 0;
		}
		else
		{
			l = -tri_nabe[3*(ltri-1)+ledg-1];
			rtri = l / 3;
			redg = (l % 3) + 1;
		}

		vbedg ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1], point_num, 
				point_xy, *tri_num, tri_vert, tri_nabe, &ltri, &ledg, &rtri, &redg );

		n = *tri_num + 1;
		l = -tri_nabe[3*(ltri-1)+ledg-1];

		for ( ; ; )
		{
			t = l / 3;
			e = ( l % 3 ) + 1;
			l = -tri_nabe[3*(t-1)+e-1];
			m2 = tri_vert[3*(t-1)+e-1];

			if ( e <= 2 )
			{
				m1 = tri_vert[3*(t-1)+e];
			}
			else
			{
				m1 = tri_vert[3*(t-1)+0];
			}

			*tri_num = *tri_num + 1;
			tri_nabe[3*(t-1)+e-1] = *tri_num;
			tri_vert[3*(*tri_num-1)+0] = m1;
			tri_vert[3*(*tri_num-1)+1] = m2;
			tri_vert[3*(*tri_num-1)+2] = m;
			tri_nabe[3*(*tri_num-1)+0] = t;
			tri_nabe[3*(*tri_num-1)+1] = *tri_num - 1;
			tri_nabe[3*(*tri_num-1)+2] = *tri_num + 1;
			top = top + 1;

			if ( point_num < top )
			{
				cout << "\n";
				cout << "DTRIS2 - Fatal error!\n";
				cout << "  Stack overflow.\n";
				delete [] stack;
				return 8;
			}

			stack[top-1] = *tri_num;

			if ( t == rtri && e == redg )
			{
				break;
			}

		}

		tri_nabe[3*(ltri-1)+ledg-1] = -3 * n - 1;
		tri_nabe[3*(n-1)+1] = -3 * (*tri_num) - 2;
		tri_nabe[3*(*tri_num-1)+2] = -l;
		ltri = n;
		ledg = 2;

		error = swapec ( m, &top, &ltri, &ledg, point_num, point_xy, *tri_num, 
				tri_vert, tri_nabe, stack );

		if ( error != 0 )
		{
			cout << "\n";
			cout << "DTRIS2 - Fatal error!\n";
			cout << "  Error return from SWAPEC.\n";
			delete [] stack;
			return error;
		}

	}
	//
	//  Now account for the sorting that we did.
	//
	for ( i = 0; i < 3; i++ )
	{
		for ( j = 0; j < *tri_num; j++ )
		{
			tri_vert[i+j*3] = indx [ tri_vert[i+j*3] - 1 ];
		}
	}

	perm_inverse ( point_num, indx );

	r82vec_permute ( point_num, indx, base, point_xy );

	delete [] indx;
	delete [] stack;

	return 0;
}

//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
	int value;

	if ( i2 < i1 )
	{
		value = i1;
	}
	else
	{
		value = i2;
	}
	return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
	int value;

	if ( i1 < i2 )
	{
		value = i1;
	}
	else
	{
		value = i2;
	}
	return value;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
	int value;

	if ( j == 0 )
	{
		cout << "\n";
		cout << "I4_MODP - Fatal error!\n";
		cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
		exit ( 1 );
	}

	value = i % j;

	if ( value < 0 )
	{
		value = value + abs ( j );
	}

	return value;
}
//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
	int value;

	if ( i < 0 )
	{
		value = -1;
	}
	else
	{
		value = 1;
	}
	return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
	int jhi;
	int jlo;
	int value;
	int wide;

	jlo = i4_min ( ilo, ihi );
	jhi = i4_max ( ilo, ihi );

	wide = jhi + 1 - jlo;

	if ( wide == 1 )
	{
		value = jlo;
	}
	else
	{
		value = jlo + i4_modp ( ival - jlo, wide );
	}

	return value;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
		int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

	int i;
	int i2hi;
	int i2lo;
	int j;
	int j2hi;
	int j2lo;

	cout << "\n";
	cout << title << "\n";
	//
	//  Print the columns of the matrix, in strips of INCX.
	//
	for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
	{
		i2hi = i2lo + INCX - 1;
		i2hi = i4_min ( i2hi, m );
		i2hi = i4_min ( i2hi, ihi );

		cout << "\n";
		//
		//  For each row I in the current range...
		//
		//  Write the header.
		//
		cout << "  Row: ";
		for ( i = i2lo; i <= i2hi; i++ )
		{
			cout << setw(6) << i << "  ";
		}
		cout << "\n";
		cout << "  Col\n";
		cout << "\n";
		//
		//  Determine the range of the rows in this strip.
		//
		j2lo = i4_max ( jlo, 1 );
		j2hi = i4_min ( jhi, n );

		for ( j = j2lo; j <= j2hi; j++ )
		{
			//
			//  Print out (up to INCX) entries in column J, that lie in the current strip.
			//
			cout << setw(5) << j << "  ";
			for ( i = i2lo; i <= i2hi; i++ )
			{
				cout << setw(6) << a[i-1+(j-1)*m] << "  ";
			}
			cout << "\n";
		}
	}

	return;
# undef INCX
}


//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
	int i;
	int value;

	if ( n <= 0 )
	{
		return 0;
	}

	value = a[0];

	for ( i = 1; i < n; i++ )
	{
		if ( a[i] < value )
		{
			value = a[i];
		}
	}
	return value; 
}
//****************************************************************************80

int lrline ( double xu, double yu, double xv1, double yv1, double xv2, 
		double yv2, double dv )

//****************************************************************************80
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
	double dx;
	double dxu;
	double dy;
	double dyu;
	double t;
	double tol;
	double tolabs;
	int value;

	tol = 100.0 * r8_epsilon ( );

	dx = xv2 - xv1;
	dy = yv2 - yv1;
	dxu = xu - xv1;
	dyu = yu - yv1;

	tolabs = tol * r8_max ( fabs ( dx ), 
			r8_max ( fabs ( dy ), 
					r8_max ( fabs ( dxu ), 
							r8_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

	t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

	if ( tolabs < t )
	{
		value = 1;
	}
	else if ( -tolabs <= t )
	{
		value = 0;
	}
	else if ( t < -tolabs )
	{
		value = -1;
	}

	return value;
}
//****************************************************************************80

bool perm_check ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
	bool found;
	int i;
	int seek;

	for ( seek = base; seek < base + n; seek++ )
	{
		found = false;

		for ( i = 0; i < n; i++ )
		{
			if ( p[i] == seek )
			{
				found = true;
				break;
			}
		}

		if ( !found )
		{
			return false;
		}

	}

	return true;
}
//****************************************************************************80

void perm_inverse ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Discussion:
//
//    This algorithm assumes that the entries in the permutation vector are
//    strictly positive.  In particular, the value 0 must no occur.
//
//    When necessary, this function shifts the data temporarily so that
//    this requirement is satisfied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
	int base;
	int i;
	int i0;
	int i1;
	int i2;
	int is;
	int p_min;

	if ( n <= 0 )
	{
		cout << "\n";
		cout << "PERM_INVERSE - Fatal error!\n";
		cout << "  Input value of N = " << n << "\n";
		exit ( 1 );
	}
	//
	//  Find the least value, and shift data so it begins at 1.
	//
	p_min = i4vec_min ( n, p );
	base = 1;

	for ( i = 0; i < n; i++ )
	{
		p[i] = p[i] - p_min + base;
	}
	//
	//  Now we can safely check the permutation.
	//
	if ( !perm_check ( n, p, base ) )
	{
		cerr << "\n";
		cerr << "PERM_INVERSE - Fatal error!\n";
		cerr << "  PERM_CHECK rejects this permutation.\n";
		exit ( 1 );
	}
	//
	//  Now we can invert the permutation.
	//
	is = 1;

	for ( i = 1; i <= n; i++ )
	{
		i1 = p[i-1];

		while ( i < i1 )
		{
			i2 = p[i1-1];
			p[i1-1] = -i2;
			i1 = i2;
		}

		is = - i4_sign ( p[i-1] );
		p[i-1] = i4_sign ( is ) * abs ( p[i-1] );
	}

	for ( i = 1; i <= n; i++ )
	{
		i1 = - p[i-1];

		if ( 0 <= i1 )
		{
			i0 = i;

			for ( ; ; )
			{
				i2 = p[i1-1];
				p[i1-1] = i0;

				if ( i2 < 0 )
				{
					break;
				}
				i0 = i1;
				i1 = i2;
			}
		}
	}
	//
	//  Now we can restore the permutation.
	//
	for ( i = 0; i < n; i++ )
	{
		p[i] = p[i] + p_min - base;
	}

	return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
	double value;

	if ( 0.0 <= x )
	{
		value = x;
	} 
	else
	{
		value = - x;
	}
	return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the 
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
	double value;

	value = 1.0;

	while ( 1.0 < ( double ) ( 1.0 + value )  )
	{
		value = value / 2.0;
	}

	value = 2.0 * value;

	return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
	double value;

	value = 1.0E+30;

	return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
	double value;

	if ( y < x )
	{
		value = x;
	} 
	else
	{
		value = y;
	}
	return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
	double value;

	if ( y < x )
	{
		value = y;
	} 
	else
	{
		value = x;
	}
	return value;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
	int value;

	if ( x < 0.0 )
	{
		value = - ( int ) ( r8_abs ( x ) + 0.5 );
	}
	else
	{
		value =   ( int ) ( r8_abs ( x ) + 0.5 );
	}

	return value;
}
//****************************************************************************80

void r82vec_permute ( int n, int p[], int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PERMUTE permutes an R82VEC in place.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array. 
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
//
//    Input/output, double A[2*N], the array to be permuted.
//
{
	double a_temp[2];
	int i;
	int iget;
	int iput;
	int istart;

	if ( !perm_check ( n, p, base ) )
	{
		cerr << "\n";
		cerr << "R82VEC_PERMUTE - Fatal error!\n";
		cerr << "  PERM_CHECK rejects this permutation.\n";
		exit ( 1 );
	}
	//
	//  In order for the sign negation trick to work, we need to assume that the
	//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
	//  So temporarily add 1-BASE to each entry to force positivity.
	//
	for ( i = 0; i < n; i++ )
	{
		p[i] = p[i] + 1 - base;
	}
	//
	//  Search for the next element of the permutation that has not been used.
	//
	for ( istart = 1; istart <= n; istart++ )
	{
		if ( p[istart-1] < 0 )
		{
			continue;
		}
		else if ( p[istart-1] == istart )
		{
			p[istart-1] = - p[istart-1];
			continue;
		}
		else
		{
			a_temp[0] = a[0+(istart-1)*2];
			a_temp[1] = a[1+(istart-1)*2];
			iget = istart;
			//
			//  Copy the new value into the vacated entry.
			//
			for ( ; ; )
			{
				iput = iget;
				iget = p[iget-1];

				p[iput-1] = - p[iput-1];

				if ( iget < 1 || n < iget )
				{
					cout << "\n";
					cout << "R82VEC_PERMUTE - Fatal error!\n";
					cout << "  Entry IPUT = " << iput << " of the permutation has\n";
					cout << "  an illegal value IGET = " << iget << ".\n";
					exit ( 1 );
				}

				if ( iget == istart )
				{
					a[0+(iput-1)*2] = a_temp[0];
					a[1+(iput-1)*2] = a_temp[1];
					break;
				}
				a[0+(iput-1)*2] = a[0+(iget-1)*2];
				a[1+(iput-1)*2] = a[1+(iget-1)*2];
			}
		}
	}
	//
	//  Restore the signs of the entries.
	//
	for ( i = 0; i < n; i++ )
	{
		p[i] = - p[i];
	}
	//
	//  Restore the base of the entries.
	//
	for ( i = 0; i < n; i++ )
	{
		p[i] = p[i] - 1 + base;
	}
	return;
}
//****************************************************************************80

int *r82vec_sort_heap_index_a ( int n, int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(*,indx(*))
//
//    or explicitly, by the call
//
//      r82vec_permute ( n, indx, base, a )
//
//    after which a(*,*) is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing, 
//    1 for 1-based indexing.
//
//    Input, double A[2*N], an array to be index-sorted.
//
//    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,R8VEC_SORT_HEAP_INDEX_A(I)).
//
{
	double aval[2];
	int i;
	int *indx;
	int indxt;
	int ir;
	int j;
	int l;

	if ( n < 1 )
	{
		return NULL;
	}

	indx = new int[n];

	for ( i = 0; i < n; i++ )
	{
		indx[i] = i;
	}

	if ( n == 1 )
	{
		indx[0] = indx[0] + base;
		return indx;
	}

	l = n / 2 + 1;
	ir = n;

	for ( ; ; )
	{
		if ( 1 < l )
		{
			l = l - 1;
			indxt = indx[l-1];
			aval[0] = a[0+indxt*2];
			aval[1] = a[1+indxt*2];
		}
		else
		{
			indxt = indx[ir-1];
			aval[0] = a[0+indxt*2];
			aval[1] = a[1+indxt*2];
			indx[ir-1] = indx[0];
			ir = ir - 1;

			if ( ir == 1 )
			{
				indx[0] = indxt;
				break;
			}
		}
		i = l;
		j = l + l;

		while ( j <= ir )
		{
			if ( j < ir )
			{
				if (   a[0+indx[j-1]*2] <  a[0+indx[j]*2] ||
						( a[0+indx[j-1]*2] == a[0+indx[j]*2] &&
								a[1+indx[j-1]*2] <  a[1+indx[j]*2] ) )
				{
					j = j + 1;
				}
			}

			if (   aval[0] <  a[0+indx[j-1]*2] ||
					( aval[0] == a[0+indx[j-1]*2] &&
							aval[1] <  a[1+indx[j-1]*2] ) )
			{
				indx[i-1] = indx[j-1];
				i = j;
				j = j + j;
			}
			else
			{
				j = ir + 1;
			}
		}
		indx[i-1] = indxt;
	}
	//
	//  Take care of the base.
	//
	for ( i = 0; i < n; i++ )
	{
		indx[i] = indx[i] + base;
	}

	return indx;
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R8MAT_DATA_READ[M*N], the table data.
//
{
	bool error;
	ifstream input;
	int i;
	int j;
	string line;
	double *table;
	double *x;

	input.open ( input_filename.c_str ( ) );

	if ( !input )
	{
		cerr << "\n";
		cerr << "R8MAT_DATA_READ - Fatal error!\n";
		cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
		return NULL;
	}

	table = new double[m*n];

	x = new double[m];

	j = 0;

	while ( j < n )
	{
		getline ( input, line );

		if ( input.eof ( ) )
		{
			break;
		}

		if ( line[0] == '#' || s_len_trim ( line ) == 0 )
		{
			continue;
		}

		error = s_to_r8vec ( line, m, x );

		if ( error )
		{
			continue;
		}

		for ( i = 0; i < m; i++ )
		{
			table[i+j*m] = x[i];
		}
		j = j + 1;

	}

	input.close ( );

	delete [] x;

	return table;
}

//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
		int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

	int i;
	int i2;
	int i2hi;
	int i2lo;
	int inc;
	int j;
	int j2hi;
	int j2lo;

	cout << "\n";
	cout << title << "\n";

	for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
	{
		i2hi = i2lo + INCX - 1;
		i2hi = i4_min ( i2hi, m );
		i2hi = i4_min ( i2hi, ihi );

		inc = i2hi + 1 - i2lo;

		cout << "\n";
		cout << "  Row: ";
		for ( i = i2lo; i <= i2hi; i++ )
		{
			cout << setw(7) << i << "       ";
		}
		cout << "\n";
		cout << "  Col\n";
		cout << "\n";

		j2lo = i4_max ( jlo, 1 );
		j2hi = i4_min ( jhi, n );

		for ( j = j2lo; j <= j2hi; j++ )
		{
			cout << setw(5) << j << " ";
			for ( i2 = 1; i2 <= inc; i2++ )
			{
				i = i2lo - 1 + i2;
				cout << setw(14) << a[(i-1)+(j-1)*m];
			}
			cout << "\n";
		}
	}

	return;
# undef INCX
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
	int n;

	n = s.length ( );

	while ( 0 < n ) 
	{
		if ( s[n-1] != ' ' )
		{
			return n;
		}
		n = n - 1;
	}

	return n;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
	char c;
	int ihave;
	int isgn;
	int iterm;
	int jbot;
	int jsgn;
	int jtop;
	int nchar;
	int ndig;
	double r;
	double rbot;
	double rexp;
	double rtop;
	char TAB = 9;

	nchar = s_len_trim ( s );
	*error = false;
	r = 0.0;
	*lchar = -1;
	isgn = 1;
	rtop = 0.0;
	rbot = 1.0;
	jsgn = 1;
	jtop = 0;
	jbot = 1;
	ihave = 1;
	iterm = 0;

	for ( ; ; )
	{
		c = s[*lchar+1];
		*lchar = *lchar + 1;
		//
		//  Blank or TAB character.
		//
		if ( c == ' ' || c == TAB )
		{
			if ( ihave == 2 )
			{
			}
			else if ( ihave == 6 || ihave == 7 )
			{
				iterm = 1;
			}
			else if ( 1 < ihave )
			{
				ihave = 11;
			}
		}
		//
		//  Comma.
		//
		else if ( c == ',' || c == ';' )
		{
			if ( ihave != 1 )
			{
				iterm = 1;
				ihave = 12;
				*lchar = *lchar + 1;
			}
		}
		//
		//  Minus sign.
		//
		else if ( c == '-' )
		{
			if ( ihave == 1 )
			{
				ihave = 2;
				isgn = -1;
			}
			else if ( ihave == 6 )
			{
				ihave = 7;
				jsgn = -1;
			}
			else
			{
				iterm = 1;
			}
		}
		//
		//  Plus sign.
		//
		else if ( c == '+' )
		{
			if ( ihave == 1 )
			{
				ihave = 2;
			}
			else if ( ihave == 6 )
			{
				ihave = 7;
			}
			else
			{
				iterm = 1;
			}
		}
		//
		//  Decimal point.
		//
		else if ( c == '.' )
		{
			if ( ihave < 4 )
			{
				ihave = 4;
			}
			else if ( 6 <= ihave && ihave <= 8 )
			{
				ihave = 9;
			}
			else
			{
				iterm = 1;
			}
		}
		//
		//  Exponent marker.
		//
		else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
		{
			if ( ihave < 6 )
			{
				ihave = 6;
			}
			else
			{
				iterm = 1;
			}
		}
		//
		//  Digit.
		//
		else if ( ihave < 11 && '0' <= c && c <= '9' )
		{
			if ( ihave <= 2 )
			{
				ihave = 3;
			}
			else if ( ihave == 4 )
			{
				ihave = 5;
			}
			else if ( ihave == 6 || ihave == 7 )
			{
				ihave = 8;
			}
			else if ( ihave == 9 )
			{
				ihave = 10;
			}

			ndig = ch_to_digit ( c );

			if ( ihave == 3 )
			{
				rtop = 10.0 * rtop + ( double ) ndig;
			}
			else if ( ihave == 5 )
			{
				rtop = 10.0 * rtop + ( double ) ndig;
				rbot = 10.0 * rbot;
			}
			else if ( ihave == 8 )
			{
				jtop = 10 * jtop + ndig;
			}
			else if ( ihave == 10 )
			{
				jtop = 10 * jtop + ndig;
				jbot = 10 * jbot;
			}

		}
		//
		//  Anything else is regarded as a terminator.
		//
		else
		{
			iterm = 1;
		}
		//
		//  If we haven't seen a terminator, and we haven't examined the
		//  entire string, go get the next character.
		//
		if ( iterm == 1 || nchar <= *lchar + 1 )
		{
			break;
		}

	}
	//
	//  If we haven't seen a terminator, and we have examined the
	//  entire string, then we're done, and LCHAR is equal to NCHAR.
	//
	if ( iterm != 1 && (*lchar) + 1 == nchar )
	{
		*lchar = nchar;
	}
	//
	//  Number seems to have terminated.  Have we got a legal number?
	//  Not if we terminated in states 1, 2, 6 or 7!
	//
	if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
	{
		*error = true;
		return r;
	}
	//
	//  Number seems OK.  Form it.
	//
	if ( jtop == 0 )
	{
		rexp = 1.0;
	}
	else
	{
		if ( jbot == 1 )
		{
			rexp = pow ( 10.0, jsgn * jtop );
		}
		else
		{
			rexp = jsgn * jtop;
			rexp = rexp / jbot;
			rexp = pow ( 10.0, rexp );
		}

	}

	r = isgn * rexp * rtop / rbot;

	return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
	int begin;
	bool error;
	int i;
	int lchar;
	int length;

	begin = 0;
	length = s.length ( );
	error = 0;

	for ( i = 0; i < n; i++ )
	{
		rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

		if ( error )
		{
			return error;
		}
		begin = begin + lchar;
		length = length - lchar;
	}

	return error;
}

//****************************************************************************80

int swapec ( int i, int *top, int *btri, int *bedg, int point_num, 
		double point_xy[], int tri_num, int tri_vert[], int tri_nabe[], 
		int stack[] )

//****************************************************************************80
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int *TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int *BTRI, *BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the points.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input/output, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//    May be updated on output because of swaps.
//
//    Input/output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list;
//    negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
	int a;
	int b;
	int c;
	int e;
	int ee;
	int em1;
	int ep1;
	int f;
	int fm1;
	int fp1;
	int l;
	int r;
	int s;
	int swap;
	int t;
	int tt;
	int u;
	double x;
	double y;
	//
	//  Determine whether triangles in stack are Delaunay, and swap
	//  diagonal edge of convex quadrilateral if not.
	//
	x = point_xy[2*(i-1)+0];
	y = point_xy[2*(i-1)+1];

	for ( ; ; )
	{
		if ( *top <= 0 ) 
		{
			break;
		}

		t = stack[(*top)-1];
		*top = *top - 1;

		if ( tri_vert[3*(t-1)+0] == i )
		{
			e = 2;
			b = tri_vert[3*(t-1)+2];
		}
		else if ( tri_vert[3*(t-1)+1] == i )
		{
			e = 3;
			b = tri_vert[3*(t-1)+0];
		}
		else
		{
			e = 1;
			b = tri_vert[3*(t-1)+1];
		}

		a = tri_vert[3*(t-1)+e-1];
		u = tri_nabe[3*(t-1)+e-1];

		if ( tri_nabe[3*(u-1)+0] == t )
		{
			f = 1;
			c = tri_vert[3*(u-1)+2];
		}
		else if ( tri_nabe[3*(u-1)+1] == t )
		{
			f = 2;
			c = tri_vert[3*(u-1)+0];
		}
		else
		{
			f = 3;
			c = tri_vert[3*(u-1)+1];
		}

		swap = diaedg ( x, y, 
				point_xy[2*(a-1)+0], point_xy[2*(a-1)+1],
				point_xy[2*(c-1)+0], point_xy[2*(c-1)+1],
				point_xy[2*(b-1)+0], point_xy[2*(b-1)+1] );

		if ( swap == 1 )
		{
			em1 = i4_wrap ( e - 1, 1, 3 );
			ep1 = i4_wrap ( e + 1, 1, 3 );
			fm1 = i4_wrap ( f - 1, 1, 3 );
			fp1 = i4_wrap ( f + 1, 1, 3 );

			tri_vert[3*(t-1)+ep1-1] = c;
			tri_vert[3*(u-1)+fp1-1] = i;
			r = tri_nabe[3*(t-1)+ep1-1];
			s = tri_nabe[3*(u-1)+fp1-1];
			tri_nabe[3*(t-1)+ep1-1] = u;
			tri_nabe[3*(u-1)+fp1-1] = t;
			tri_nabe[3*(t-1)+e-1] = s;
			tri_nabe[3*(u-1)+f-1] = r;

			if ( 0 < tri_nabe[3*(u-1)+fm1-1] )
			{
				*top = *top + 1;
				stack[(*top)-1] = u;
			}

			if ( 0 < s )
			{
				if ( tri_nabe[3*(s-1)+0] == u )
				{
					tri_nabe[3*(s-1)+0] = t;
				}
				else if ( tri_nabe[3*(s-1)+1] == u )
				{
					tri_nabe[3*(s-1)+1] = t;
				}
				else
				{
					tri_nabe[3*(s-1)+2] = t;
				}

				*top = *top + 1;

				if ( point_num < *top )
				{
					return 8;
				}

				stack[(*top)-1] = t;
			}
			else
			{
				if ( u == *btri && fp1 == *bedg )
				{
					*btri = t;
					*bedg = e;
				}

				l = - ( 3 * t + e - 1 );
				tt = t;
				ee = em1;

				while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
				{
					tt = tri_nabe[3*(tt-1)+ee-1];

					if ( tri_vert[3*(tt-1)+0] == a )
					{
						ee = 3;
					}
					else if ( tri_vert[3*(tt-1)+1] == a )
					{
						ee = 1;
					}
					else
					{
						ee = 2;
					}
				}
				tri_nabe[3*(tt-1)+ee-1] = l;
			}

			if ( 0 < r )
			{
				if ( tri_nabe[3*(r-1)+0] == t )
				{
					tri_nabe[3*(r-1)+0] = u;
				}
				else if ( tri_nabe[3*(r-1)+1] == t )
				{
					tri_nabe[3*(r-1)+1] = u;
				}
				else
				{
					tri_nabe[3*(r-1)+2] = u;
				}
			}
			else
			{
				if ( t == *btri && ep1 == *bedg )
				{
					*btri = u;
					*bedg = f;
				}

				l = - ( 3 * u + f - 1 );
				tt = u;
				ee = fm1;

				while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
				{
					tt = tri_nabe[3*(tt-1)+ee-1];

					if ( tri_vert[3*(tt-1)+0] == b )
					{
						ee = 3;
					}
					else if ( tri_vert[3*(tt-1)+1] == b )
					{
						ee = 1;
					}
					else
					{
						ee = 2;
					}

				}
				tri_nabe[3*(tt-1)+ee-1] = l;
			}
		}
	}

	return 0;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 29

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time ( NULL );
	tm = localtime ( &now );

	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	if ( len != 0 )
	{
		cout << time_buffer << "\n";
	}

	return;
# undef TIME_SIZE
}
//****************************************************************************80

void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num, 
		int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg )

//****************************************************************************80
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 December 2008
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Modified:
//
//    02 September 2003
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//
//    Input, int TRI_NABE[TRI_NUM*3], the triangle neighbor list; negative
//    values are used for links of a counter clockwise linked list of boundary
//    edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int *RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int *REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
	int a;
	double ax;
	double ay;
	int b;
	double bx;
	double by;
	bool done;
	int e;
	int l;
	int lr;
	int t;
	//
	//  Find the rightmost visible boundary edge using links, then possibly
	//  leftmost visible boundary edge using triangle neighbor information.
	//
	if ( *ltri == 0 )
	{
		done = false;
		*ltri = *rtri;
		*ledg = *redg;
	}
	else
	{
		done = true;
	}

	for ( ; ; )
	{
		l = -tri_nabe[3*((*rtri)-1)+(*redg)-1];
		t = l / 3;
		e = 1 + l % 3;
		a = tri_vert[3*(t-1)+e-1];

		if ( e <= 2 )
		{
			b = tri_vert[3*(t-1)+e];
		}
		else
		{
			b = tri_vert[3*(t-1)+0];
		}

		ax = point_xy[2*(a-1)+0];
		ay = point_xy[2*(a-1)+1];

		bx = point_xy[2*(b-1)+0];
		by = point_xy[2*(b-1)+1];

		lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

		if ( lr <= 0 )
		{
			break;
		}

		*rtri = t;
		*redg = e;

	}

	if ( done )
	{
		return;
	}

	t = *ltri;
	e = *ledg;

	for ( ; ; )
	{
		b = tri_vert[3*(t-1)+e-1];
		e = i4_wrap ( e-1, 1, 3 );

		while ( 0 < tri_nabe[3*(t-1)+e-1] )
		{
			t = tri_nabe[3*(t-1)+e-1];

			if ( tri_vert[3*(t-1)+0] == b )
			{
				e = 3;
			}
			else if ( tri_vert[3*(t-1)+1] == b )
			{
				e = 1;
			}
			else
			{
				e = 2;
			}

		}

		a = tri_vert[3*(t-1)+e-1];
		ax = point_xy[2*(a-1)+0];
		ay = point_xy[2*(a-1)+1];

		bx = point_xy[2*(b-1)+0];
		by = point_xy[2*(b-1)+1];

		lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

		if ( lr <= 0 )
		{
			break;
		}

	}

	*ltri = t;
	*ledg = e;

	return;
}
