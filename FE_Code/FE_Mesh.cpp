#include "FE_Mesh.h"

// Other headers
#include <iostream>
#include <algorithm>
using namespace std;

// -------------------------------------------------------------
// -------------------------------------------------------------
//   class FE_Mesh
// -------------------------------------------------------------
// -------------------------------------------------------------

FE_Mesh::FE_Mesh():
nNodes_(0)
{

}

FE_Mesh::~FE_Mesh()
{

}

// -------------------------------------------------------------
//   mappings from element id to associated nodes
// -------------------------------------------------------------
void FE_Mesh::element_connectivity_global(const int eltID,
		Array<int> & nodes) const
{
//	const int npe = get_nNodesPerElement();
	const int npe = 4;
	nodes.reset(npe);

	// use connectivity arrays
	for (int inode = 0; inode < npe; inode++)
		nodes(inode) = connectivity_(inode, eltID);
}

// -------------------------------------------------------------
//   element_coordinates
// -------------------------------------------------------------
void FE_Mesh::element_coordinates(const int eltID,
		DENS_MAT & xCoords) const
{
//	const int nne = get_nNodesPerElement();
	const int nne = 4;

	xCoords.reset(nSD_, nne, false);
	for (int inode=0; inode<nne; inode++)
	{
		const int id = connectivity_(inode, eltID);
		for (int isd=0; isd<nSD_; isd++)
			xCoords(isd,inode) = nodalCoords_(isd,id);
	}
}


// -------------------------------------------------------------
// Nodes on a given plane (x=xv) */
// -------------------------------------------------------------
int FE_Mesh::get_Node_xy(double xv, double yv)
{

	int inode;
	for (inode=0; inode<nNodes_; inode++)
	{
		if ( ( std::fabs(nodalCoords_(0, inode)- xv)<std::numeric_limits<double>::epsilon() ) &&
		     ( std::fabs(nodalCoords_(1, inode)- yv)<std::numeric_limits<double>::epsilon() ) )
			break;
	}
	return(inode);
}


// -------------------------------------------------------------
// Nodes on a given plane (x=xv) */
// -------------------------------------------------------------
void FE_Mesh::get_Nodes_x(double xv, DENS_VEC & nodes) const
{

	int nNds = 0;
	for (int inode=0; inode<nNodes_; inode++)
	{
		if ( std::fabs(nodalCoords_(0, inode)- xv)<std::numeric_limits<double>::epsilon() ) nNds++;
	}

	nodes.reset(nNds);
	nNds = 0;
	for (int inode=0; inode<nNodes_; inode++)
	{
		if ( std::fabs(nodalCoords_(0, inode)- xv)<std::numeric_limits<double>::epsilon() )  nodes(nNds++) = inode;
	}

}

// -------------------------------------------------------------
// Nodes on a given plane (y=yv) */
// -------------------------------------------------------------
void FE_Mesh::get_Nodes_y(double yv, DENS_VEC & nodes) const
{
	int nNds = 0;
	for (int inode=0; inode<nNodes_; inode++)
	{
		if ( std::fabs(nodalCoords_(1, inode)- yv)<std::numeric_limits<double>::epsilon() ) nNds++;
	}

	nodes.reset(nNds);
	nNds = 0;
	for (int inode=0; inode<nNodes_; inode++)
	{
		if ( std::fabs(nodalCoords_(1, inode)- yv)<std::numeric_limits<double>::epsilon() ) nodes(nNds++) = inode;
	}
}


// -------------------------------------------------------------
// Nodes on a given area
// -------------------------------------------------------------
void FE_Mesh::get_Nodes(DENS_VEC indpd, DENS_VEC coords, double * xpd, DENS_VEC &nodespd, DENS_VEC &nodes)
{

	int nNds = 0;

	for (int inode=0; inode < nds_overlp_.size(); inode++)
	{
		for (int ielem=0; ielem<indpd.size();ielem++)
		{
			if ( ( std::fabs(nodalCoords_(0, nds_overlp_(inode))- coords(2*ielem+0))<1e-9 ) &&
			     ( std::fabs(nodalCoords_(1, nds_overlp_(inode))- coords(2*ielem+1))<1e-9 ) ){
				nNds++;
			}
		}
	}
	nodes.reset(nNds);
	nodespd.reset(nNds);

	nNds = 0;

	// garanties that node match
	for (int ielem=0; ielem<indpd.size();ielem++)
	{
		for (int inode=0; inode < nds_overlp_.size(); inode++)
		{
			if ( ( std::fabs(nodalCoords_(0, nds_overlp_(inode))- coords(2*ielem+0))<1e-9 ) &&
				 ( std::fabs(nodalCoords_(1, nds_overlp_(inode))- coords(2*ielem+1))<1e-9 ) ){
				nodes(nNds) = nds_overlp_(inode);
				nodespd(nNds++) = indpd(ielem);

				int offset = 3*indpd(ielem);
				xpd[offset+0] = nodalCoords_(0, nds_overlp_(inode));
				xpd[offset+1] = nodalCoords_(1, nds_overlp_(inode));
			}
		}
	}
}

// -------------------------------------------------------------
// Interpolate nodes on a given element
// -------------------------------------------------------------
void FE_Mesh::compute_shape_functions()
{
	// local coords of nodes: localCoords_(isd, ip)
	DENS_MAT localCoords_;
	localCoords_.resize(2,4);

//	// Matrix of local nodal coordinates
//	localCoords_(0,0) = -1; localCoords_(1,0) = -1;
//	localCoords_(0,1) = +1; localCoords_(1,1) = -1;
//	localCoords_(0,2) = +1; localCoords_(1,2) = +1;
//	localCoords_(0,3) = -1; localCoords_(1,3) = +1;

	// Compute shape functions at ip's
	// Set sizes of matrices and vectors
	if ((int)N_.size() != nElts_) N_.resize(nElts_);
	for (int isd = 0; isd < nElts_; isd++) N_[isd].resize(max_nPDNds_Elt_, 4);

	// Compute shape functions at interpolated points
	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt) {

		// Matrix of local nodal coordinates
		localCoords_(0,0) = localx_(0, ielt); localCoords_(1,0) = localy_(0, ielt);
		localCoords_(0,1) = localx_(1, ielt); localCoords_(1,1) = localy_(1, ielt);
		localCoords_(0,2) = localx_(2, ielt); localCoords_(1,2) = localy_(2, ielt);
		localCoords_(0,3) = localx_(3, ielt); localCoords_(1,3) = localy_(3, ielt);

		for (int ip = 0; ip < max_nPDNds_Elt_; ip++) {
			double r = localx_(ip, ielt);
			double s = localy_(ip, ielt);

			for (int iNode = 0; iNode < 4; iNode++) {
				double rI = localCoords_(0, iNode);
				double sI = localCoords_(1, iNode);
				N_[ielt](ip, iNode) = 0.25 * (1 + r*rI) * (1 + s*sI);
			}
		}
	}

//	// Compute shape functions at interpolated points
//	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt) {
//
//		for (int ip = 0; ip < max_nPDNds_Elt_; ip++) {
//			double r = localx_(ip, ielt);
//			double s = localy_(ip, ielt);
//
//			for (int iNode = 0; iNode < 4; iNode++) {
//				double rI = localCoords_(0, iNode);
//				double sI = localCoords_(1, iNode);
//				N_[ielt](ip, iNode) = 0.25 * (1 + r*rI) * (1 + s*sI);
//			}
//		}
//	}

//	cout << " " << endl;
//	cout << " localx_.print(); " << endl;
//	localx_.print();
//
//	cout << " " << endl;
//	cout << " localy_.print(); " << endl;
//	localy_.print();
}

//{
//	// local coords of nodes: localCoords_(isd, ip)
//	DENS_MAT localCoords_;
//	localCoords_.resize(2,4);
//
////	// Matrix of local nodal coordinates
////	localCoords_(0,0) = -1; localCoords_(1,0) = -1;
////	localCoords_(0,1) = +1; localCoords_(1,1) = -1;
////	localCoords_(0,2) = +1; localCoords_(1,2) = +1;
////	localCoords_(0,3) = -1; localCoords_(1,3) = +1;
//
//	// Compute shape functions at ip's
//	// Set sizes of matrices and vectors
//	if ((int)N_.size() != nElts_) N_.resize(nElts_);
//	for (int isd = 0; isd < nElts_; isd++) N_[isd].resize(max_nPDNds_Elt_, 4);
//
//	// Compute shape functions at interpolated points
//	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt) {
//
//		// Matrix of local nodal coordinates
//		localCoords_(0,0) = localx_(0, ielt); localCoords_(1,0) = localy_(0, ielt);
//		localCoords_(0,1) = localx_(1, ielt); localCoords_(1,1) = localy_(1, ielt);
//		localCoords_(0,2) = localx_(2, ielt); localCoords_(1,2) = localy_(2, ielt);
//		localCoords_(0,3) = localx_(3, ielt); localCoords_(1,3) = localy_(3, ielt);
//
//		for (int ip = 0; ip < max_nPDNds_Elt_; ip++) {
//			double r = localx_(ip, ielt);
//			double s = localy_(ip, ielt);
//
//			for (int iNode = 0; iNode < 4; iNode++) {
//				double rI = localCoords_(0, iNode);
//				double sI = localCoords_(1, iNode);
//				N_[ielt](ip, iNode) = 0.25 * (1 + r*rI) * (1 + s*sI);
//			}
//		}
//	}
//}

void FE_Mesh::interpolate_solution(DENS_VEC pd_nds, DENS_VEC sol, double *ItrpSol)
{

	int ielt, ip1, ip3, ip, offsetpd;

	double dmpx, dmpy;
	dmpx=0; dmpy=0;

	// Compute shape functions at interpolated points
//	for (ielt = (nElts_-nElts_Overlp_+68+60); ielt < (nElts_); ++ielt) {
	for (ielt = (nElts_-nElts_Overlp_+nElts_Overlp_ext); ielt < nElts_; ++ielt) {
		for (ip1 = 0; ip1 < max_nPDNds_Elt_; ip1++) {

			offsetpd = element_points_(ip1, ielt);

//			// Pass to PD the nodes that match only
//			for (ip3 = 0; ip3 < pd_nds.size(); ip3++) {
//
//				if ( pd_nds(ip3) == offsetpd ){

					// Interpolate the values
					for (ip = 0; ip < 4; ip++) {

						const int offsetfe = connectivity_(ip, ielt);

						dmpx += sol(2*offsetfe+0)*N_[ielt](ip1, ip);
						dmpy += sol(2*offsetfe+1)*N_[ielt](ip1, ip);

					}

					ItrpSol[3*offsetpd+0] = dmpx;
					ItrpSol[3*offsetpd+1] = dmpy;
					dmpx=0; dmpy=0;
//				}
//			}
		}
	}
}


// -------------------------------------------------------------
// Find nodes on a given element
// -------------------------------------------------------------
void FE_Mesh::find_Nodes(DENS_VEC nds_overlp_pd, double * x)
{
	int nNds, i, offset, offset3d, offset2d, id, inode;
	double disv = 1e-9; // value that increases the element

	// Find maximum number of nodes in the elements of the overlapping zone
	max_nPDNds_Elt_ = 0;
//	for (int ielt=(nElts_-nElts_Overlp_); ielt<nElts_;ielt++)
	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt)
	{
		nNds = 0;
		for (int j=0; j<nds_overlp_pd.size();j++)
		{
			i = nds_overlp_pd(j);
			offset = 3*i;
			if ( ( x[offset+0]>(minx(ielt)-disv) ) && ( x[offset+1]>(miny(ielt)-disv) ) &&
				 ( x[offset+0]<(maxx(ielt)+disv) ) && ( x[offset+1]<(maxy(ielt)+disv) ) )
			{
				nNds++;
			}
		}
		max_nPDNds_Elt_ = max(max_nPDNds_Elt_, nNds);
	}

	// Allocates memory
	element_points_.reset(max_nPDNds_Elt_, nElts_);
	localx_.reset(max_nPDNds_Elt_, nElts_);
	localy_.reset(max_nPDNds_Elt_, nElts_);

	int flag = 1;
//	for (int ielt=(nElts_-nElts_Overlp_); ielt<nElts_;ielt++)
	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt)
	{
		nNds = 4;
		for (int j=0; j<nds_overlp_pd.size();j++)
		{
			flag = 1;
			i = nds_overlp_pd(j);
			offset = 3*i;
			if ( ( x[offset+0]>(minx(ielt)-disv) ) && ( x[offset+1]>(miny(ielt)-disv) ) &&
			     ( x[offset+0]<(maxx(ielt)+disv) ) && ( x[offset+1]<(maxy(ielt)+disv) ) ) {

				// This if garantees that the node matches with the element node
				// Was used to match the 11 grid mesh
				for (inode=0; inode<4;inode++){
					id = connectivity_(inode, ielt);

					// This if garantees that the node matches with the element node
					if ( (  std::fabs( x[offset+0]-nodalCoords_(0, id) ) < 1e-7)  &&  ( std::fabs( x[offset+1]-nodalCoords_(1, id) ) < 1e-7 ) ) {
						element_points_(inode, ielt) = i;
						flag = 0;
					}
				}
				if (flag==1) {
					element_points_(nNds++, ielt) = i;
				}

			}
		}
	}

	double xc,yc;
//	for (int ielt=(nElts_-nElts_Overlp_); ielt<nElts_;ielt++) {
//	for (int ielt=0; ielt<nElts_Overlp_;ielt++) {
	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt) {

		xc = (nodalCoords_(0, connectivity_(0, ielt)) + nodalCoords_(0, connectivity_(1, ielt)) + nodalCoords_(0, connectivity_(2, ielt)) + nodalCoords_(0, connectivity_(3, ielt)))/4.0;
		yc = (nodalCoords_(1, connectivity_(0, ielt)) + nodalCoords_(1, connectivity_(1, ielt)) + nodalCoords_(1, connectivity_(2, ielt)) + nodalCoords_(1, connectivity_(3, ielt)))/4.0;

		for (int i = 0; i < max_nPDNds_Elt_; ++i) {
			offset2d = 2*ielt;
			offset3d = 3*element_points_(i, ielt);

			localx_(i, ielt) = 2 * (x[offset3d+0]-xc) / std::fabs(minx(ielt)-maxx(ielt));
			localy_(i, ielt) = 2 * (x[offset3d+1]-yc) / std::fabs(miny(ielt)-maxy(ielt));

			// Takes off noise
			if (std::fabs( localx_(i, ielt) ) < 1e-7) localx_(i, ielt) = 0.0;
			if (std::fabs( localy_(i, ielt) ) < 1e-7) localy_(i, ielt) = 0.0;

		}
	}

//	cout << " " << endl;
//	cout << " elements points " << endl;
//	element_points_.print();
//
//	cout << " " << endl;
//	cout << " localx_ " << endl;
//	localx_.print();
//
//	cout << " " << endl;
//	cout << " localy_ " << endl;
//	localy_.print();

	compute_shape_functions();

}

// -------------------------------------------------------------
// -------------------------------------------------------------
//   class FE_Uniform2DMesh
// -------------------------------------------------------------
// -------------------------------------------------------------

FE_Uniform2DMesh::FE_Uniform2DMesh(FILE* File)
{
	int NmbrElements, NmbrNodes;
	int e1, e2, e3, e4;
	double ax, ay;
	char buf[100];

	nx_[0] = 11;
	nx_[1] = 11;

	borders_[0][0] = -1.0;
	borders_[1][0] =  1.0;
	borders_[0][1] = -1.0;
	borders_[1][1] =  1.0;
//	cout << " " << endl;
//	cout << " maxx" << endl;
//	maxx.print();
//
//	cout << " " << endl;
//	cout << " maxy" << endl;
//	maxy.print();

	xscale_ = 1.0;
	yscale_ = 1.0;
	zscale_ = 1.0;

	// Member data setup
	nSD_ = 2;

	fgets(buf,100,File);

	if (strstr(buf,"COORDINATES") ){
		//Reading the number of nodes.
		fscanf(File, "%i", &NmbrNodes);
		nNodes_ = NmbrNodes;
		nodalCoords_.reset(2, NmbrNodes);

		// Fill nodal coordinates
		for (int inode = 0; inode < NmbrNodes; ++inode) {
			//	Reading the number of nodes and coordinates.
			fscanf(File, "%lf %lf", &ax, &ay);
			//	Assigns 'x' and 'y' coordinates.
			nodalCoords_(0,inode) = ax;
			nodalCoords_(1,inode) = ay;
		}
	}
//	nodalCoords_.print();

	fscanf(File, "%lf", &ax);
	fgets(buf,100,File);
	if(strstr(buf,"CONNECTIVITY")) {
		//Reading the number of elements.
		fscanf(File, "%i", &NmbrElements);
		nElts_ = NmbrElements;
		connectivity_.reset(4, NmbrElements);

		// Fill element connectivities
		for (int ielt = 0; ielt < NmbrElements; ++ielt) {
			fscanf(File, "%i %i %i %i", &e1, &e2, &e3, &e4);
			connectivity_(0,ielt) = e1-1;
			connectivity_(1,ielt) = e2-1;
			connectivity_(2,ielt) = e3-1;
			connectivity_(3,ielt) = e4-1;
		}
	}
//	connectivity_.print();
	fscanf(File, "%lf", &ax);
	fgets(buf,100,File);
	if(strstr(buf,"OVERLAPPING")) {
		fscanf(File, "%i", &e1);
		nElts_Overlp_ = e1;
	}
//	nodalCoords_.print();

	//	connectivity_.print();
	fscanf(File, "%lf", &ax);
	fgets(buf,100,File);
	if(strstr(buf,"MINUS")) {
		fscanf(File, "%i", &e1);
		nElts_Overlp_ext = e1-1;
	}
	//	nodalCoords_.print();
	cout << " NUMBER MINUS NODES ELEMENTS " << nElts_Overlp_ext << endl;


	int nNds = 0;
	int iflag = 0;
	DENS_VEC nds_overlp_aux;

	nds_overlp_aux.reset(4*nElts_Overlp_);

//	for (int ielt=(nElts_-nElts_Overlp_); ielt<nElts_;ielt++)
	for (int ielt = (nElts_-nElts_Overlp_); ielt < nElts_; ++ielt)
	{
		for (int ienode=0; ienode<4;ienode++)
		{
			const int id = connectivity_(ienode, ielt);
			nds_overlp_aux(nNds++) = id+1;
		}
	}
	//nds_overlp_aux.print();

	nNds = 0;
	iflag = 0;
	for (int ielt=0; ielt<4*nElts_Overlp_;ielt++)
	{
		for (int ienode=0; ienode<ielt;ienode++)
		{
			iflag = 0;
			if (nds_overlp_aux(ielt) == nds_overlp_aux(ienode))
			{
				iflag = 1;
				break;
			}
		}
		if (iflag == 0)	nNds++;
	}

	nds_overlp_.reset(nNds);
	nNds = 0;
	iflag = 0;
	for (int ielt=0; ielt<4*nElts_Overlp_;ielt++)
	{
		for (int ienode=0; ienode<ielt;ienode++)
		{
			iflag = 0;
			if (nds_overlp_aux(ielt) == nds_overlp_aux(ienode))
			{
				iflag = 1;
				break;
			}
		}
		if (iflag == 0)	nds_overlp_(nNds++) = nds_overlp_aux(ielt)-1;
	}

	/** Min and Max of the elements in the overlapping */
	minx.reset(nElts_);
	maxx.reset(nElts_);
	miny.reset(nElts_);
	maxy.reset(nElts_);

	double id1, id2, idi;

	// find min and max
	for (int ielt=0; ielt<nElts_;ielt++)
	{
		// min x
		id1 = nodalCoords_(0, connectivity_(0, ielt));
		id2 = nodalCoords_(0, connectivity_(1, ielt));
		idi = min(id1,id2);
		id1 = nodalCoords_(0, connectivity_(2, ielt));
		id2 = min(id1,idi);
		id1 = nodalCoords_(0, connectivity_(3, ielt));
		idi = min(id1,id2);
		minx(ielt) = idi;

		// min y
		id1 = nodalCoords_(1, connectivity_(0, ielt));
		id2 = nodalCoords_(1, connectivity_(1, ielt));
		idi = min(id1,id2);
		id1 = nodalCoords_(1, connectivity_(2, ielt));
		id2 = min(id1,idi);
		id1 = nodalCoords_(1, connectivity_(3, ielt));
		idi = min(id1,id2);
		miny(ielt) = idi;

		// max x
		id1 = nodalCoords_(0, connectivity_(0, ielt));
		id2 = nodalCoords_(0, connectivity_(1, ielt));
		idi = max(id1,id2);
		id1 = nodalCoords_(0, connectivity_(2, ielt));
		id2 = max(id1,idi);
		id1 = nodalCoords_(0, connectivity_(3, ielt));
		idi = max(id1,id2);
		maxx(ielt) = idi;

		// max y
		id1 = nodalCoords_(1, connectivity_(0, ielt));
		id2 = nodalCoords_(1, connectivity_(1, ielt));
		idi = max(id1,id2);
		id1 = nodalCoords_(1, connectivity_(2, ielt));
		id2 = max(id1,idi);
		id1 = nodalCoords_(1, connectivity_(3, ielt));
		idi = max(id1,id2);
		maxy(ielt) = idi;
	}

	/** Alpha and Beta for the interpolation F = Alpha*t + Beta */
	alph_bt_x.reset(2*nElts_);
	alph_bt_y.reset(2*nElts_);

	localx_.reset(4, NmbrElements);
	localy_.reset(4, NmbrElements);


	/** Min and Max of the elements in the overlapping */
	DENS_VEC minxx;
	DENS_VEC maxxx;
	DENS_VEC minyx;
	DENS_VEC maxyx;
	minxx.reset(nElts_);
	maxxx.reset(nElts_);
	minyx.reset(nElts_);
	maxyx.reset(nElts_);


	// find min and max
	for (int ielt=0; ielt<nElts_;ielt++)
	{
		// min x
		id1 = nodalCoords_(0, connectivity_(0, ielt));
		minxx(ielt) = id1;

		// min y
		id1 = nodalCoords_(1, connectivity_(0, ielt));
		minyx(ielt) = id1;

		// max x
		id1 = nodalCoords_(0, connectivity_(2, ielt));
		maxxx(ielt) = id1;

		// max y
		id1 = nodalCoords_(1, connectivity_(2, ielt));
		maxyx(ielt) = id1;
	}

	for (int ielt = 0; ielt < nElts_; ++ielt) {

		double mmx = maxxx(ielt)-minxx(ielt);
		alph_bt_x(2*ielt+0) = 2.0/mmx;
		alph_bt_x(2*ielt+1) = -(maxxx(ielt)+minxx(ielt))/mmx;

		double mmy = maxyx(ielt)-minyx(ielt);
		alph_bt_y(2*ielt+0) = 2.0/mmy;
		alph_bt_y(2*ielt+1) = -(maxyx(ielt)+minyx(ielt))/mmy;
	}

//	cout << " " << endl;
//	cout << " alph_bt_x.print(); " << endl;
//	alph_bt_x.print();
//
//	cout << " " << endl;
//	cout << " alph_bt_y.print(); " << endl;
//	alph_bt_y.print();
//	cout << " " << endl;

}

FE_Uniform2DMesh::FE_Uniform2DMesh(const int nx,
									const int ny,
									double xmin, double xmax,
									double ymin, double ymax,
									double xscale,
									double yscale)
{
	nx_[0] = nx;
	nx_[1] = ny;

	borders_[0][0] = xmin;
	borders_[1][0] = xmax;
	borders_[0][1] = ymin;
	borders_[1][1] = ymax;

	xscale_ = xscale;
	yscale_ = yscale;

	// Compute region size and element size
	for (int i = 0; i < 2; i++) {
		Lx_[i] = borders_[1][i] - borders_[0][i];
		dx_[i] = Lx_[i]/nx_[i];
	}

	// Member data setup
	nSD_    = 2;
	nElts_  = nx_[0] * nx_[1];
	nNodes_ = (nx_[0]+1) * (nx_[1]+1);

	nodalCoords_.reset(2, nNodes_);
    connectivity_.reset(4, nElts_);

	// Fill nodal coordinates
	double ix[2];
	int inode = 0;
	for (int j = 0; j <= nx_[1]; ++j) {
		ix[1] = borders_[0][1] + j*dx_[1];
		for (int i = 0; i <= nx_[0]; ++i) {
			ix[0] = borders_[0][0] + i*dx_[0];
			for (int m = 0; m < 2; ++m) {
				nodalCoords_(m,inode) = ix[m];
			}
			++inode;
		}
	}

	// Compute element connectivities
	int ielt = 0;
	int noffx = 1;
	int noffy = nx_[0] + 1;
	for (int j = 0; j < nx_[1]; ++j) {
		for (int i = 0; i < nx_[0]; ++i) {
			int i1 = i + j*noffy;
			connectivity_(0,ielt) = i1;
			connectivity_(1,ielt) = i1 + noffx;
			connectivity_(2,ielt) = i1 + noffx + noffy;
			connectivity_(3,ielt) = i1 + noffy;
			++ielt;
		}
	}
}

FE_Uniform2DMesh::~FE_Uniform2DMesh()
{
	// Nothing to do here
}
