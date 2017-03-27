#ifndef FE_ELEMENT_H
#define FE_ELEMENT_H

#include "MatrixLibrary.h"
#include "Array2D.h"

#include <vector>

/**
 *  @class  FE_Element
 *  @brief  Base class for a finite element holding info for canonical element
 */
class FE_Element {

public:
	FE_Element();

	FE_Element(int nSD,
			   int nDofs,
			   int numEltNodes,
			   int numIPs,
			   int numFaces,
			   int numFaceNodes,
			   int numFaceIPs);

	virtual ~FE_Element();

    /** get number of element nodes */
    int num_elt_nodes() { return numEltNodes_; }

    /** get number of integration points */
    int num_ips() { return numIPs_; }

    /** get number of degrees of freedom */
    int num_node_dofs() {return 2;}

    DENS_MAT shape_functions() {return N_;}

    vector<DENS_MAT> d_shape_functions() {return dN_;}

	virtual void set_quadrature(int quadrature_type) = 0;

	enum FE_ElementQuadrature { GAUSSIAN_QUADRATURE, NODAL_QUADRATURE};



protected:

	// Number of spatial dimensions
	int nSD_;

	// Number of degree of freedom
	int nDofs_;

	// Number of element nodes
	int numEltNodes_;

	// Number of integration points
	int numIPs_;

	// local coords of nodes: localCoords_(isd, ip)
	DENS_MAT localCoords_;

	// quadrature scheme: localCoords_(isd, ip)
	DENS_MAT ipCoords_; // local coordinates

	// matrix of shape functions at ip's: N_(ip, node)
	DENS_MAT N_;

	vector<DENS_MAT> dN_;

	// integration point weights: ipWeights_(ip)
	DENS_VEC ipWeights_; // local coordinates

};


/**
 *  @class  FE_ElementQuad
 *  @author Fabiano F. Bargos Jan 25/2012
 *  @brief  2D, linear 4-node quadrilateral element with 2x2 quadrature
 */
class FE_ElementQuad : public FE_Element {

public:
	FE_ElementQuad();
	~FE_ElementQuad();

	void element_stiffness(double **** dsde, DENS_MAT &coord, DENS_MAT &el_mat);

protected:
	virtual void set_quadrature(int quadrature_type);
};

#endif // FE_ELEMENT_H
