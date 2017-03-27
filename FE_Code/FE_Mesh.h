#ifndef FE_MESH_H
#define FE_MESH_H

#include "Array.h"
#include "Array2D.h"
#include "MatrixLibrary.h"

// Other headers
#include <vector>
#include <map>
#include <set>
#include <utility>

#include <cmath>
#include <limits>

using namespace std;

/**
 *  @class  FE_Mesh
 *  @brief  Base class for a finite element mesh
 */
class FE_Mesh {

public:
	/** constructor */
	FE_Mesh();

	/** destructor */
	~FE_Mesh();

	/** return connectivity (global ID numbers) for element eltID */
	void element_connectivity_global(const int eltID,
			Array<int> & nodes) const;

	/**
	 *  return spatial coordinates for element nodes on eltID,
	 *  indexed xCoords(isd,inode)
	 */
	void element_coordinates(const int eltID,
			DENS_MAT & xCoords) const;

	/** access to the nodal coordinate values */
	const DENS_MAT & nodal_coordinates(void) {return nodalCoords_  ;}

	/** access to the element connectivity values */
	const Array2D<int> & connectivity(void) {return connectivity_  ;}
		
	/** return number of spatial dimensions */
	int get_nSpatialDimensions() const { return nSD_; };

	/** return total number of nodes */
	int get_nNodes() const { return nNodes_; };

	/** return number of elements */
	int get_nElements() const { return nElts_; };

	/** return the nodes on a desired x value */
	void get_Nodes_x(double xv, DENS_VEC & nodes) const;

	/** return the nodes on a desired y value */
	void get_Nodes_y(double yv, DENS_VEC & nodes) const;

	/** return the nodes on a desired x value */
	int get_Node_xy(double xv, double yv);

	/** return the nodes on a desired overlapping area */
	void get_Nodes(DENS_VEC indpd, DENS_VEC coords, double * x, DENS_VEC &nodespd, DENS_VEC &nodes);

	void find_Nodes(DENS_VEC nds_overlp_pd, double * x);

	void compute_shape_functions();

	void interpolate_solution(DENS_VEC nds_interp_pd, DENS_VEC sol, double *ItrpSol);

protected:

	/** number of spatial dimensions */
	int nSD_;

	/** number of elements */
	int nElts_;

	/** number of elements in the overlapping area*/
	int nElts_Overlp_;

	/** number of nodes */
	int nNodes_;

	/** Nodal coordinates: nodalCoords_(nsd, numnode) */
	DENS_MAT nodalCoords_;

	/** Element connectivity: connectivity_(neltnode, nelt) */
	Array2D<int> connectivity_;

	/** length scaling used by lammps */
	double xscale_, yscale_, zscale_;

	/** Nodes in the overlapping area */
	DENS_VEC nds_overlp_;

	/** Min and Max of the elements in the overlapping */
	DENS_VEC minx;
	DENS_VEC maxx;
	DENS_VEC miny;
	DENS_VEC maxy;

	/** Alpha and Beta for the interpolation F = Alpha*t + Beta */
	DENS_VEC alph_bt_x;
	DENS_VEC alph_bt_y;

	/** local coordinates local_x_(neltnode, nelt) */
	Array2D<double> localx_;
	Array2D<double> localy_;

	/** Shape functions computed on the collocation points */
	vector<DENS_MAT> N_;

	/** Points inside the element: element_points_(neltnode, nelt) */
	Array2D<int> element_points_;

	/** Number Points inside the element */
	int max_nPDNds_Elt_;

	/** number of elements in the overlapping area (external part)*/
	int nElts_Overlp_ext;
};

/**
 *  @class  FE_Uniform2DMesh
 *  @author Fabiano F. Bargos Jan 25/2012(int nNodes, DENS_VEC indpd, DENS_VEC coords, DENS_VEC &nodespd, DENS_VEC &nodes)
 *  @brief  Derived class for a uniform mesh, a structured mesh with
 *          fixed element sizes in x, y, and z directions.
 */

class FE_Uniform2DMesh : public FE_Mesh {

public:
	/** constructor */
	FE_Uniform2DMesh(FILE* File);

	/** constructor */
	FE_Uniform2DMesh(const int nx,
					 const int ny,
					 double xmin, double xmax,
					 double ymin, double ymax,
					 double xscale,
					 double yscale);

	/** destructor */
	~FE_Uniform2DMesh();


protected:

	// Number of elements in each spatial direction
	int nx_[2];

	// Bounds of region on which mesh is defined
	double borders_[2][2];

	// Region size in each direction
	double Lx_[2];

	// Element size in each direction
	double dx_[2];

	/** create global-to-unique node mapping */
	void setup_periodicity();

};


#endif // FE_MESH_H
