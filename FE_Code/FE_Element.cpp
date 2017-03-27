#include "FE_Element.h"

// Other headers
#include "math.h"

#include <iostream>
#include <cstdlib>
using namespace std;

// -------------------------------------------------------------
// -------------------------------------------------------------
//   class FE_Element
// -------------------------------------------------------------
// -------------------------------------------------------------
FE_Element::FE_Element()
{

}

FE_Element::FE_Element(int nSD,
		               int nDofs,
					   int numEltNodes,
					   int numIPs,
					   int numFaces,
					   int numFaceNodes,
					   int numFaceIPs)
: nSD_(nSD),
  nDofs_(nDofs),
  numEltNodes_(numEltNodes),
  numIPs_(numIPs),
  localCoords_(nSD,numEltNodes),
  ipCoords_(nSD,numIPs),
  N_(numIPs,numEltNodes),
  dN_(nSD),
  ipWeights_(numIPs)
  {

  }

FE_Element::~FE_Element()
{
	// Nothing else to do
}


// -------------------------------------------------------------
// -------------------------------------------------------------
//   class FE_ElementQuad
// -------------------------------------------------------------
// -------------------------------------------------------------

FE_ElementQuad::FE_ElementQuad()
: FE_Element(2, 2, 4, 4, 1, 2, 2)
{

	// Matrix of local nodal coordinates
	localCoords_(0,0) = -1; localCoords_(1,0) = -1;
	localCoords_(0,1) = +1; localCoords_(1,1) = -1;
	localCoords_(0,2) = +1; localCoords_(1,2) = +1;
	localCoords_(0,3) = -1; localCoords_(1,3) = +1;

	//       3 ----- 2     y
	//      /       /      |
	//     /       /       |
	//    0 ----- 1         ---> x
	//                    /
	//                   /
	//                  z

	set_quadrature(GAUSSIAN_QUADRATURE);

}

//-----------------------------------------------------------------
void FE_ElementQuad::set_quadrature(int quadrature)
{
	double a = 1./sqrt(3.);
	if (quadrature == NODAL_QUADRATURE) {
		a = 1.0;
	}
	// Matrix of integration point locations  (Gaussian) & follows local conn
	ipCoords_(0,0) = -a; ipCoords_(1,0) = -a;
	ipCoords_(0,1) = +a; ipCoords_(1,1) = -a;
	ipCoords_(0,2) = +a; ipCoords_(1,2) = +a;
	ipCoords_(0,3) = -a; ipCoords_(1,3) = +a;


	// Compute shape functions at ip's
	for (int ip = 0; ip < numIPs_; ip++) {
		double r = ipCoords_(0, ip);
		double s = ipCoords_(1, ip);

		for (int iNode = 0; iNode < numEltNodes_; iNode++) {
			double rI = localCoords_(0, iNode);
			double sI = localCoords_(1, iNode);

			N_(ip, iNode) = 0.25 * (1 + r*rI) * (1 + s*sI);
		}
	}

	// Compute shape functions at ip's
	// Set sizes of matrices and vectors
	if ((int)dN_.size() != nSD_) dN_.resize(nSD_);
	for (int isd = 0; isd < nSD_; isd++) dN_[isd].resize(numIPs_, numEltNodes_);

	// Loop over integration points
	for (int ip = 0; ip < numIPs_; ip++) {
		double r = ipCoords_(0, ip);
		double s = ipCoords_(1, ip);

		for (int iNode = 0; iNode < numEltNodes_; iNode++) {
			double rI = localCoords_(0, iNode);
			double sI = localCoords_(1, iNode);

			dN_[0](ip, iNode) = 0.25 * rI * (1 + s*sI);
			dN_[1](ip, iNode) = 0.25 * sI * (1 + r*rI);
		}
	}

	// Integration point weights
	ipWeights_ = 1.0;
	
}

//DENS_MAT coords;
//element_coordinates(eltID, &coords);
//-----------------------------------------------------------------
void FE_ElementQuad::element_stiffness(double **** dsde, DENS_MAT &coord, DENS_MAT &el_mat)
//
//  Assemble the element stiffness
//
//    Arguments;
//
//      nSD_               # coordinates (2 or 3 for 2D or 3D problem)
//      nDofs_             # degrees of freedom per node (not often nDofs_ =  nSD_)
//      numEltNodes_       # nodes on the element
//      coords(i,a)        ith coord of ath node
//
//   Local variables
//      numIPs_            # integration points
//      w(intpt)           weight for integration point no. intpt
//      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
//      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
//      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
//      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
//      det                Determinant of jacobian
//      dsde[i][j][k][l]   Derivative of stress_ij with respect:strain_kl
//      el_mat(row,col)       Rows && cols of element stiffness
//
{
	DENS_MAT dNdx(nSD_, numEltNodes_);
	DENS_MAT dxdxi(nSD_, nSD_);
	//DENS_MAT el_mat(nDofs_* numEltNodes_, nDofs_*numEltNodes_);

	//
	//  Set up integration points && weights
	//
	//xilist = integrationumIPs_( nSD_,numEltNodes_,numIPs_, elident);
	//w = integrationweights( nSD_, numEltNodes_,numIPs_,elident);

	DENS_MAT dNdxi(nDofs_, numEltNodes_);
	DENS_MAT dxidx(nDofs_, numEltNodes_);

	//
	//  Loop over the integration points
	//
	for(int intpt = 0; intpt < numIPs_; ++intpt) {
		//
		//     Compute derivatives wrt local coords
		//
		// copy derivatives on int points
		for(int j = 0; j <  nSD_; ++j) {
			for(int a = 0; a < numEltNodes_; ++a) {
				dNdxi(j, a) = dN_[j](intpt, a);
			}
		}

		//
		//     Compute the jacobian matrix && its determinant
		//
		for(int i = 0; i <  nSD_; ++i) {
			for(int j = 0; j <  nSD_; ++j) {
				dxdxi(i,j) = 0.;
				for(int a = 0; a < numEltNodes_; ++a) {
					dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(j, a);
				}
			}
		}

		//DENS_MAT dxidx = inv(dxdxi);
		double dt  = ( dxdxi(0,0)*dxdxi(1,1) - dxdxi(1,0)*dxdxi(0,1) );
		dxidx(0,0) =  dxdxi(1, 1);
		dxidx(0,1) = -dxdxi(0, 1);
		dxidx(1,0) = -dxdxi(1, 0);
		dxidx(1,1) =  dxdxi(0, 0);
		dxidx      = (1.0/dt)*dxidx;

		//
		//     Convert shape function derivatives:derivatives wrt global coords
		//
		for(int a = 0; a < numEltNodes_; ++a) {
			for(int i = 0; i <  nSD_; ++i) {
				dNdx(i, a) = 0.;
				for(int j = 0; j <  nSD_; ++j) {
					dNdx(i, a) = dNdx(i, a) + dNdxi(j, a) * dxidx(i, j);
				}
			}
		}
		//
		//     Compute the element stiffness
		//
		int row, col;
		for(int a = 0; a < numEltNodes_; ++a) {
			for(int i = 0; i < nDofs_; ++i) {
				for(int b = 0; b < numEltNodes_; ++b) {
					for(int k = 0; k < nDofs_; ++k) {
						row = nDofs_*(a)+i;
						col = nDofs_*(b)+k;
						for(int j = 0; j <  nSD_; ++j) {
							for(int l = 0; l <  nSD_; ++l) {
								//cout << " row " << " " << row << " col " << " " << col << endl;
//								el_mat(col,row) = el_mat(col,row) + dsde[i][j][k][l]*dNdx(b,l)*dNdx(a,j)*ipWeights_(intpt)*dt;
								el_mat(col,row) = el_mat(col,row) + dsde[i][j][k][l]*dNdx(l, b)*dNdx(j, a)*1*dt;
							}
						}
					}
				}
			}
		}
	}
}



FE_ElementQuad::~FE_ElementQuad()
{
	// Nothing to do here
}
