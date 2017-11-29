/* ----------------------------------------------------------------------
   Drive program for coupling Peridynamics/Finite Element
   Fabiano F. Bargos
   Brown University/State University of Campinas
   November/2011

   This routine links LAMMPS
   (Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories)
   to a generic Finite Element Code.

   The implementation is based on the "c++_driver" available in LAMMPS.

   c++_driver = simple example of how an umbrella program
   can invoke LAMMPS as a library on some subset of procs
Syntax: c++_driver P in.lammps
P = # of procs to run LAMMPS on
must be <= # of procs the driver code itself runs on
in.lammps = LAMMPS input script
See README for compilation instructions
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

#include "FE_Engine.h"
#include "MatrixLibrary.h"

#include <iostream>
using namespace std;
using namespace LAMMPS_NS;

#ifndef BCVALUE
#define BCVALUE 0.56
#endif

#ifndef TLOAD
#define TLOAD 1
#endif

// -------------------------------------------------------------
// Define functions */
// -------------------------------------------------------------
void ReadFEandPDFiles(double *PDCoords, DENS_VEC &lp_ext_nds_bc2PD,
    DENS_VEC &lp_int_nds_bc2FE,
    DENS_VEC &lp_ext_coo_bc2PD,
    DENS_VEC &lp_int_coo_bc2FE,
    DENS_VEC &lp_ov_nds);

// -------------------------------------------------------------
// Main Program */
// -------------------------------------------------------------
int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);

  // -------------------------------------------------------------
  // Lammps solver */
  // -------------------------------------------------------------
  int me = 0;
  int lammps = 1;

  // open LAMMPS input script

  FILE *fp;
  if (me == 0) {
    //		fp = fopen(arg[1],"r");
    fp = stdin;
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line

  // Treating input args to pass to LAMMPS
  int narg_copy = 1;
  char** arg_copy;
  arg_copy    = (char**)malloc(1*sizeof(char*));
  arg_copy[0] =  (char*)malloc(100*sizeof(char));
  strcpy(arg_copy[0],arg[0]);


  //	LAMMPS *lmp = new LAMMPS(narg, arg, MPI_COMM_WORLD);
  LAMMPS *lmp = new LAMMPS(narg_copy, arg_copy, MPI_COMM_WORLD);

  if (lammps == 1)  lmp->input->file();

  int natoms = static_cast<int> (lmp->atom->natoms);

  // Stores the initial position of the atoms
  double *x0 = new double[3*natoms];
  lammps_get_coords(lmp,x0);          // no LAMMPS class function for this

  // this is for convergence
  double *xn   = new double[3*natoms];

  // external loop convergence
  double *xen = new double[3*natoms];
  double *xe0 = new double[3*natoms];


  // Stores the position of the atoms
  double *x = new double[3*natoms];
  lammps_get_coords(lmp,x);          // no LAMMPS class function for this

  cout << "LAMMPS was set" << endl;
  fclose (fp);

  // -------------------------------------------------------------
  // Finite element solver */
  // -------------------------------------------------------------

  // Material definition
  int ndof   = 2;
  int ncoord = 2;

  // Should be used for non-damaged problems
  double mu  = 10;
  double nu  = 0.33333333333333333333333;
  // Load
  double load = 0.2e0/5.;

//  double mu   = 72e9;
//  double nu   = 0.33333333333333333333333;
////  double load = 2.0e8; // no damage
//  double load = 2.5e7;
//  double topload = 5e7;

//  double analA = load/0.05/mu*(1.-nu)*(1.+nu);
//  double analB = -nu*analA/(1.-nu);


  int flag = atoi(arg[2]); // uses 1 dirichet; 2 for neumann; 3 robin
  int thetaflag = 0; //uses 0 without aitkenl 1 for aitken

  FILE *File;

  // Open test file
  if((File = fopen(arg[1], "r")) == NULL) {
    cout << "Error in opening" << arg[0] << endl;
  }

  ofstream lammpsfile;
  lammpsfile.open ("lammpsconv.txt");
  ofstream errorfile;
  errorfile.open ("errorconv.txt");

  FE_Engine feEngine(File);
  feEngine.FEsetup(mu, nu, ndof, ncoord, load);

  fclose (File);
  cout << "FE was set" << endl;

  /** CREATES A FE MESH ON PD REGION (FABIANO JUNE 27/2016) */
  int number_dim = 2;
  DENS_MAT coords_xy0;
  DENS_MAT PDConnectivity;
  coords_xy0.reset(number_dim, natoms);

  int offset3d = 0;
  for(int i = 0; i < natoms; ++i) {
    coords_xy0(0, i) = x0[offset3d + 0];
    coords_xy0(1, i) = x0[offset3d + 1];
    offset3d += 3;
  }

  int number_triangles;
  feEngine.GenerateMesh(coords_xy0, PDConnectivity);

  // -------------------------------------------------------------
  // Reads information of the overlapping */
  // -------------------------------------------------------------
  DENS_VEC lp_extrnl_nds_bc2PD;		// Stores the BC nodes
  DENS_VEC lp_extrnl_nds_bc2PDnew;	// Stores the BC nodes
  DENS_VEC lp_intrnl_nds_bc2FE;		// Stores the BC nodes
  DENS_VEC lp_intrnl_nds_bc2FEnew;	// Stores the BC nodes
  DENS_VEC lp_extrnl_coo_bc2PD;		// Stores the BC nodes coordinates
  DENS_VEC lp_intrnl_coo_bc2FE;		// Stores the BC nodes coordinates
  DENS_VEC fe_nds_bc2PD;			// Stores the BC nodes
  DENS_VEC fe_nds_bc2FE;			// Stores the BC nodes
  DENS_VEC fe_nds_bc2FEinv;                     // Stores the BC nodes
  DENS_VEC dirichletpresc_values;		// Stores the BC nodes
  DENS_VEC neumann_values;		// Stores the BC nodes
  DENS_VEC lp_overl_nds;			// Stores the nodes in the overlapping

  ReadFEandPDFiles(x, lp_extrnl_nds_bc2PD,
      lp_intrnl_nds_bc2FE,
      lp_extrnl_coo_bc2PD,
      lp_intrnl_coo_bc2FE,
      lp_overl_nds);

  // -------------------------------------------------------------
  // Matching FE and PD overlapping area */
  // -------------------------------------------------------------
  fprintf ( stdout, "\n\n Finding FE nodes in the BC2PD region... \n" );
  feEngine.feMesh->get_Nodes(lp_extrnl_nds_bc2PD, lp_extrnl_coo_bc2PD, x, lp_extrnl_nds_bc2PDnew, fe_nds_bc2PD);

  fprintf ( stdout, "\n\n Finding FE nodes in the BC2FE region... \n" );
  feEngine.feMesh->get_Nodes(lp_intrnl_nds_bc2FE, lp_intrnl_coo_bc2FE, x, lp_intrnl_nds_bc2FEnew, fe_nds_bc2FE);

  // Set for each element which are the peridynamics nodes inside
  feEngine.feMesh->find_Nodes(lp_overl_nds, x);

  // -------------------------------------------------------------
  // Initial FE/PD loop */
  // -------------------------------------------------------------
  int offset2d, idrch = 0, idrch2 = 0;
  DENS_VEC upCoords, upCoords2, fesol, fesol2;

  double * QkFE  = new double [fe_nds_bc2FE.size()*2];
  double * Qk1FE = new double [fe_nds_bc2FE.size()*2];
  for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
    QkFE[i*2]   = 0.0;
    QkFE[i*2+1] = 0.0;
  }
  double Qkin, Qksq;

  // -------------------------------------------------------------
  // used to compute neumann bc */
  // -------------------------------------------------------------
  neumann_values.reset(ndof * fe_nds_bc2FE.size());

  // -------------------------------------------------------------
  // used to compute dirichlet bc */
  // -------------------------------------------------------------
  dirichletpresc_values.reset(ndof * fe_nds_bc2FE.size());

  DENS_MAT Coords;
  Coords =  feEngine.feMesh->nodal_coordinates();

  double ElemEdge = abs( Coords(0, 1) - Coords(0, 0) );
  double rcoeff=atof(arg[3]); //Robin
  
  lammpsfile << "case: " << flag << " rcoef: " << rcoeff << endl;

  double analA = load/ElemEdge/mu*(1.-nu)*(1.+nu);
  double analB = -nu*analA/(1.-nu);

  fprintf(stdout, "ElemEdge=%lf\n analA=%lf, analB=%lf\n", ElemEdge, analA, analB);

  double *xdisp    = new double[3*natoms];
  double *xdisp1   = new double[3*natoms];
  double *xdispold = new double[3*natoms];
  double *QkPD     = new double[3*natoms];
  double *Qk1PD    = new double[3*natoms];
  double *QkPDall  = new double[3*natoms];
  double *Qk1PDall = new double[3*natoms];

  for (int i = 0; i < 3*natoms; i++) {
    xdisp[i]    = 0.0; // displacement from FE interpolated to PD
    xdisp1[i]   = 0.0;
    QkPD[i]     = 0.0;
    Qk1PD[i]    = 0.0;
    QkPDall[i]  = 0.0;
    Qk1PDall[i] = 0.0;
  }
  double QkinPD, QksqPD;

  int ntstp         = 10000;
  double theta      = 0.0;
  double thetapd    = 0.0;
  if(thetaflag) thetapd=0.5;
  double thetapdall = 0.0;  

  // -------------------------------------------------------------
  // Main FE/PD loop */
  // -------------------------------------------------------------
  int mytsp = 0;
  int tstp = 0;
  double CoupleTol = 1e-15;
  double diffx     = 100.0;
  double diffx0    = 100.0;
  double diffFE    = 100.0;
  double diffFE0   = 100.0;
  double analerror, lammperror;
  double initFE = 100.0;

  DENS_MAT Neighbor;
  Neighbor.reset(3, fe_nds_bc2FE.size());
  DENS_MAT Normal;
  Normal.reset(2, fe_nds_bc2FE.size());

  fe_nds_bc2FEinv.reset(fe_nds_bc2FE.size());

  // For Neumann and Robin only
  if (flag != 1) {

    for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
      int NodeID = fe_nds_bc2FE(i);
      cout << i << " vs " << fe_nds_bc2FE(i) << endl;

      fe_nds_bc2FEinv(fe_nds_bc2FE(i)) = i;

      if ( fe_nds_bc2FE(i) == 0 ) {
        Neighbor(0, NodeID) = 0;
        Neighbor(1, NodeID) = 1;
        Neighbor(2, NodeID) = fe_nds_bc2FE.size()-1;
      }
      else if ( fe_nds_bc2FE(i) == (fe_nds_bc2FE.size()-1) ) {
        Neighbor(0, NodeID) = fe_nds_bc2FE.size()-1;
        Neighbor(1, NodeID) = 0;
        Neighbor(2, NodeID) = fe_nds_bc2FE.size()-2;
      }
      else{
        Neighbor(0, NodeID) = NodeID;
        Neighbor(1, NodeID) = NodeID+1;
        Neighbor(2, NodeID) = NodeID-1;
      }
      // EDGE 1 y = 0.6
      if ( (Coords(1, NodeID)) > BCVALUE ){ 
        Normal(0, NodeID) = 0.0; Normal(1, NodeID) = -1.0;
      }
      // EDGE 2 x = 0.6
      else if  ( (Coords(0, NodeID)) > BCVALUE ) {
        Normal(0, NodeID) = -1.0; Normal(1, NodeID) = 0.0;
      }
      // EDGE 3 y = -0.6
      else if  ( (Coords(1, NodeID)) < -BCVALUE ){
        Normal(0, NodeID) = 0.0; Normal(1, NodeID) = 1.0;
      }
      // EDGE 4 x = -0.6
      else if  ( (Coords(0, NodeID)) < -BCVALUE ){
        Normal(0, NodeID) = 1.0; Normal(1, NodeID) = 0.0;
      }
    }
    for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
      int NodeID = fe_nds_bc2FE(i);
      cout << Neighbor(0,NodeID) << " " << Neighbor(1,NodeID) << " " << Neighbor(2,NodeID) << " " << Normal(0, NodeID) << " " << Normal(1, NodeID) << endl;
    }
  }

  // -------------------------------------------------------------
  // USED TO compute NEUMANN bc (FE) */
  // -------------------------------------------------------------
  DENS_MAT coords_xy;
  DENS_MAT Stress;
  coords_xy.reset(number_dim, natoms);

  lammps_get_coords(lmp,x);
  for(int i = 0; i < 3*natoms; ++i) xe0[i] = x[i];

for (int loadi=0; loadi<TLOAD; loadi++)
{
//  feEngine.FEsetup1(mu, nu, ndof, ncoord, load+(topload-load)/TLOAD*(loadi+1));
  tstp = 0;
  diffFE    = 100.0;
  diffFE0    = 100.0;
  while((tstp < ntstp)&&((diffFE/diffFE0 > CoupleTol) || (tstp < 2))) {

    diffx  = 0.0;
    diffFE = 0.0;

    // no LAMMPS class function for this
    lammps_get_coords(lmp,x);

    for(int i = 0; i < 3*natoms; ++i) xen[i] = x[i];

    // -------------------------------------------------------------
    // Compute dirichlet bc (FE) */
    // -------------------------------------------------------------
    if (flag == 1) { 		

      idrch = 0;
      for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
        offset3d = 3 * lp_intrnl_nds_bc2FEnew(i);
        offset2d = ndof * fe_nds_bc2FE(i);
        dirichletpresc_values(idrch) = (x[offset3d+0] - Coords(0, fe_nds_bc2FE(i)));
        idrch++;
        dirichletpresc_values(idrch) = (x[offset3d+1] - Coords(1, fe_nds_bc2FE(i)));
        idrch++;
      }
    }  

    // -------------------------------------------------------------
    // Compute Neumann and Robin bc (FE) */
    // -------------------------------------------------------------
    if (flag == 2) {
      // computing stress matrix - (in the rows) Sx, Sy, Sxy - the elements are in the cols
      offset3d = 0;
      for(int i = 0; i < natoms; ++i) {
        coords_xy(0, i) = x[offset3d + 0];        
        coords_xy(1, i) = x[offset3d + 1];        
        offset3d += 3;
      }	  
      feEngine.ComputeStress(nu, mu, coords_xy0, coords_xy, PDConnectivity, Stress);

      // computing neumann Rhs
      idrch = 0;
      for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {

        int PDNode = lp_intrnl_nds_bc2FEnew(i);
        int FENode = fe_nds_bc2FE(i);

        int NLt = Neighbor(1,FENode);
        int NRt = Neighbor(2,FENode);

        int PDNLt = lp_intrnl_nds_bc2FEnew(fe_nds_bc2FEinv(NLt));
        int PDNRt = lp_intrnl_nds_bc2FEnew(fe_nds_bc2FEinv(NRt));

        //int PDNLt = lp_intrnl_nds_bc2FEnew(NLt);
        //int PDNRt = lp_intrnl_nds_bc2FEnew(NRt);

        double Sx, Sy;
        Sx = ((2.0/3.0) * (Stress(0, PDNode) * Normal(0, FENode) + Stress(2, PDNode) * Normal(1, FENode)) +
            (1.0/6.0) * (Stress(0, PDNLt)  * Normal(0, FENode) + Stress(2, PDNLt)  * Normal(1, FENode)) +
            (1.0/6.0) * (Stress(0, PDNRt)  * Normal(0, FENode) + Stress(2, PDNRt)  * Normal(1, FENode)) );
        Sy = ((2.0/3.0) * (Stress(2, PDNode) * Normal(0, FENode) + Stress(1, PDNode) * Normal(1, FENode)) +
            (1.0/6.0) * (Stress(2, PDNLt)  * Normal(0, FENode) + Stress(1, PDNLt) *  Normal(1, FENode)) +
            (1.0/6.0) * (Stress(2, PDNRt)  * Normal(0, FENode) + Stress(1, PDNRt) *  Normal(1, FENode)) );

        // At the corners
//        if ( (FENode == 0) || (FENode == 24) || (FENode == 48) || (FENode == 72) )
       if ( ( fabs( fabs(Coords(0, FENode)) - fabs(Coords(1, FENode)) ) < 1.0e-8 ) &&
              (fabs(Coords(0, FENode)) > BCVALUE) && (fabs(Coords(1, FENode)) > BCVALUE) ) 
       {        
          Sx = ((1.0/3.0) * (Stress(0, PDNode) * Normal(0, NRt) + Stress(2, PDNode) * Normal(1, NRt)) +
              (1.0/6.0) * (Stress(0, PDNRt)  * Normal(0, NRt) + Stress(2, PDNRt)  * Normal(1, NRt)) +					      
              (1.0/3.0) * (Stress(0, PDNode) * Normal(0, NLt) + Stress(2, PDNode) * Normal(1, NLt)) +
              (1.0/6.0) * (Stress(0, PDNLt)  * Normal(0, NLt) + Stress(2, PDNLt)  * Normal(1, NLt)) );

          Sy = ((1.0/3.0) * (Stress(2, PDNode) * Normal(0, NRt) + Stress(1, PDNode) * Normal(1, NRt)) +
              (1.0/6.0) * (Stress(2, PDNRt)  * Normal(0, NRt) + Stress(1, PDNRt)  * Normal(1, NRt)) +					      
              (1.0/3.0) * (Stress(2, PDNode) * Normal(0, NLt) + Stress(1, PDNode) * Normal(1, NLt)) +
              (1.0/6.0) * (Stress(2, PDNLt)  * Normal(0, NLt) + Stress(1, PDNLt)  * Normal(1, NLt)) );
        }

        // x direction
        neumann_values(idrch) = ElemEdge * Sx;		
        idrch++;

        // y direction
        neumann_values(idrch) = ElemEdge * Sy;
        idrch++;
      }

      cout << " NEUMANN VALUES " << "with elmt size="<< ElemEdge <<  endl;
      idrch = 0;
      for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
        cout << neumann_values(idrch) << " ";
        idrch++;
        cout << neumann_values(idrch) << endl;
        idrch++;
      }


      idrch = 0;
      // Neumann
      if (flag == 2) {
        for(int i = 0; i < 2*fe_nds_bc2FE.size(); ++i) {
          dirichletpresc_values(idrch) = neumann_values(idrch); 
          idrch++;
        }
      }
    }
    // Robin
    else if (flag == 3){

      // computing stress matrix - (in the rows) Sx, Sy, Sxy - the elements are in the cols
      offset3d = 0;
      for(int i = 0; i < natoms; ++i) {
        coords_xy(0, i) = x[offset3d + 0];        
        coords_xy(1, i) = x[offset3d + 1];        
        offset3d += 3;
      }	  
      feEngine.ComputeStress(nu, mu, coords_xy0, coords_xy, PDConnectivity, Stress);

      // computing neumann Rhs
      idrch = 0;
      for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {

        int PDNode = lp_intrnl_nds_bc2FEnew(i);
        int FENode = fe_nds_bc2FE(i);

        int NLt = Neighbor(1,FENode);
        int NRt = Neighbor(2,FENode);

        int PDNLt = lp_intrnl_nds_bc2FEnew(fe_nds_bc2FEinv(NLt));
        int PDNRt = lp_intrnl_nds_bc2FEnew(fe_nds_bc2FEinv(NRt));

        //int PDNLt = lp_intrnl_nds_bc2FEnew(NLt);
        //int PDNRt = lp_intrnl_nds_bc2FEnew(NRt);

        offset3d = 3 * lp_intrnl_nds_bc2FEnew(i);

        double Sx, Sy;
        Sx = ((2.0/3.0) * ((Stress(0, PDNode) * Normal(0, FENode) + Stress(2, PDNode) * Normal(1, FENode)) + 
              rcoeff * (x[3*PDNode+0] - Coords(0, FENode))) +
            (1.0/6.0) * ((Stress(0, PDNLt)  * Normal(0, FENode) + Stress(2, PDNLt)  * Normal(1, FENode)) + 
              1.*rcoeff * (x[3*PDNLt+0] - Coords(0, NLt))) +
            (1.0/6.0) * ((Stress(0, PDNRt)  * Normal(0, FENode) + Stress(2, PDNRt)  * Normal(1, FENode)) +
              1.*rcoeff * (x[3*PDNRt+0] - Coords(0, NRt))));
        Sy = ((2.0/3.0) * ((Stress(2, PDNode) * Normal(0, FENode) + Stress(1, PDNode) * Normal(1, FENode)) +
              rcoeff * (x[3*PDNode+1] - Coords(1, FENode))) +
            (1.0/6.0) * ((Stress(2, PDNLt)  * Normal(0, FENode) + Stress(1, PDNLt) *  Normal(1, FENode)) +
              1.*rcoeff * (x[3*PDNLt+1] - Coords(1, NLt))) +
            (1.0/6.0) * ((Stress(2, PDNRt)  * Normal(0, FENode) + Stress(1, PDNRt) *  Normal(1, FENode)) +
              1.*rcoeff * (x[3*PDNRt+1] - Coords(1, NRt))));

        // At the corners
//        if ( (FENode == 0) || (FENode == 24) || (FENode == 48) || (FENode == 72) )
	if ( ( fabs( fabs(Coords(0, FENode)) - fabs(Coords(1, FENode)) ) < 1.0e-8 ) &&
              (fabs(Coords(0, FENode)) > BCVALUE) && (fabs(Coords(1, FENode)) > BCVALUE) )
       {
          Sx = ((1.0/3.0) * ((Stress(0, PDNode) * Normal(0, NRt) + Stress(2, PDNode) * Normal(1, NRt)) +
                rcoeff * (x[3*PDNode+0] - Coords(0, FENode))) +
              (1.0/6.0) * ((Stress(0, PDNRt)  * Normal(0, NRt) + Stress(2, PDNRt)  * Normal(1, NRt)) +	
                1.*rcoeff * (x[3*PDNRt+0] - Coords(0, NRt))) +
              (1.0/3.0) * ((Stress(0, PDNode) * Normal(0, NLt) + Stress(2, PDNode) * Normal(1, NLt)) +
                rcoeff * (x[3*PDNode+0] - Coords(0, FENode))) +
              (1.0/6.0) * ((Stress(0, PDNLt)  * Normal(0, NLt) + Stress(2, PDNLt)  * Normal(1, NLt)) +
                1.*rcoeff * (x[3*PDNLt+0] - Coords(0, NLt))));

          Sy = ((1.0/3.0) * ((Stress(2, PDNode) * Normal(0, NRt) + Stress(1, PDNode) * Normal(1, NRt)) +
                rcoeff * (x[3*PDNode+1] - Coords(1, FENode))) +
              (1.0/6.0) * ((Stress(2, PDNRt)  * Normal(0, NRt) + Stress(1, PDNRt)  * Normal(1, NRt)) +			      
                1.*rcoeff * (x[3*PDNRt+1] - Coords(1, NRt))) +
              (1.0/3.0) * ((Stress(2, PDNode) * Normal(0, NLt) + Stress(1, PDNode) * Normal(1, NLt)) +
                rcoeff * (x[3*PDNode+1] - Coords(1, FENode))) +
              (1.0/6.0) * ((Stress(2, PDNLt)  * Normal(0, NLt) + Stress(1, PDNLt)  * Normal(1, NLt)) +
                1.*rcoeff * (x[3*PDNLt+1] - Coords(1, NLt))));
        }

        // x direction
        neumann_values(idrch) = ElemEdge * Sx;		
        idrch++;

        // y direction
        neumann_values(idrch) = ElemEdge * Sy;
        idrch++;
      }

      cout << " ROBIN VALUES " << "with elmt size="<< ElemEdge <<  endl;
      idrch = 0;
      for(int i = 0; i < fe_nds_bc2FE.size(); ++i) {
        cout << neumann_values(idrch) << " ";
        idrch++;
        cout << neumann_values(idrch) << endl;
        idrch++;
      }

      idrch = 0;
      for(int i = 0; i < 2*fe_nds_bc2FE.size(); ++i) {
        dirichletpresc_values(idrch) = neumann_values(idrch); 
        idrch++;
      }
    }

    // -------------------------------------------------------------
    // Solving FE */
    // -------------------------------------------------------------
    //feEngine.FESolve(mytsp++, fe_nds_bc2FE, dirichletpresc_values, upCoords, fesol, flag);
    feEngine.FESolve(mytsp++, fe_nds_bc2FE, dirichletpresc_values, upCoords, fesol, flag, rcoeff*ElemEdge, Neighbor);


    // -------------------------------------------------------------
    // Compute dirichlet bc (Lammps) */
    // -------------------------------------------------------------
    for (int i = 0; i < 3*natoms; i++) xdisp[i] = 0.0; // displacement from FE interpolated to PD

    feEngine.feMesh->interpolate_solution(lp_extrnl_nds_bc2PD, fesol, xdisp);

    QkinPD = 0.0;
    QksqPD = 0.0;
    for(int i = 0; i < 3*natoms; ++i) {
      Qk1PD[i] = QkPD[i];
      QkPD[i]  = xdisp1[i] - xdisp[i];
      QkinPD  += (Qk1PD[i]-QkPD[i])*QkPD[i];
      QksqPD  += (Qk1PD[i]-QkPD[i])*(Qk1PD[i]-QkPD[i]);
    }

    thetapd += (thetapd-1.)*QkinPD/(QksqPD+1.e-10);

    if (thetapd<-1.) thetapd=-1.;
    if (thetapd>0.95) thetapd=0.95;

    if(!thetaflag)
      thetapd = 0.0; // dirichlet-dirichlet

    for(int i = 0; i < 3*natoms; ++i) {
      x[i]        = x0[i] + thetapd*xdisp1[i] + (1-thetapd)*xdisp[i];
      xdisp1[i]   = thetapd*xdisp1[i] + (1-thetapd)*xdisp[i];
      xdispold[i] = 0.0;
      Qk1PDall[i] = 0.0;
      QkPDall[i]  = 0.0;
    }
    int ntstpPDMax = 1000;
    int PDloopi    = 0.0;
    double PDtol   = 1e-12;
    double S       = 100.0;
    double S0      = S;
    thetapdall     = 0.0;

    // Update Lammps coordinates and run
    lammps_put_coords(lmp, x);
    for(int i = 0; i < 1; ++i) 
    {
      lmp->input->one("run 2000");
      //      lmp->input->one("minimize 1.0e-12 1.0e-12 10000 100000");
    }

    // no LAMMPS class function for this
    lammps_get_coords(lmp,xn);

    S = 0.0;
    for(int j = 0; j < 3*natoms; ++j) {
      S   += (xn[j] - x[j]) * (xn[j] - x[j]);
      x[j] = xn[j];
    }

    cout << "Iteration : " << PDloopi << " " << S/S0 << endl;

    for(int j = 0; j < 3*natoms; ++j) {
      diffx += (x[j] - xen[j]) * (x[j] - xen[j]);
    }
    if (tstp == 0) diffx0 = diffx;

    if (tstp == 0) fesol2.reset( fesol.size() );

    diffFE = 0.0;
    initFE = 0.0;
    for (int i = 0; i < fesol.size(); ++i) {
      diffFE   += (fesol2(i) - fesol(i)) * (fesol2(i) - fesol(i)) ;
      initFE   += fesol(i)*fesol(i);
      fesol2(i) = fesol(i);
    }

    analerror = 0.0;
    for (int i=0;(2*i) < fesol.size();i++)
    {
      analerror += (fesol(2*i)-analA*(Coords(0, i)+1.))*(fesol(2*i)-analA*(Coords(0, i)+1.)) + (fesol(2*i+1)-analB*(Coords(1, i)+1.))*(fesol(2*i+1)-analB*(Coords(1, i)+1.));
    }
    lammperror = 0.0;
    for(int i = 0; i < natoms; ++i)
    {
	lammperror += (x[3*i]-analA*(xe0[3*i]+1.)-xe0[3*i])*(x[3*i]-analA*(xe0[3*i]+1.)-xe0[3*i]) + (x[3*i+1]-analB*(xe0[3*i+1]+1.)-xe0[3*i+1])*(x[3*i+1]-analB*(xe0[3*i+1]+1.)-xe0[3*i+1]);
    }


    if (tstp >= 1) diffFE0 = initFE;	

    lammpsfile << "Iteration : " << tstp << " Lammps steps: " << PDloopi << " Lammps Error: " << S/S0;
    lammpsfile << " Loop Error: " << diffFE/diffFE0 << " theta: " << theta << " thetaPD: " << thetapd << endl;
//    lammpsfile << " fesolsize: " << fesol.size() <<  " fesol1: " << fesol(2*1132) << " fesol2: " << fesol(2*1052) << " afesol1: " << analA*(Coords(0, 1132)) << " afesol2: " << analA*(Coords(0, 1052)) << endl;

    errorfile << analerror << " " << lammperror << endl;

    tstp++;
  }
}

  // -------------------------------------------------------------
  // Ending codes */
  // -------------------------------------------------------------
  delete [] x;
  delete [] x0;
  delete [] xdisp;
  delete [] xn;
  delete [] xen;
  delete [] xe0;

  if (lammps == 1) delete lmp;

  printf("\n Ending test... \n");

  lammpsfile.close();
  MPI_Finalize();
}


// -------------------------------------------------------------
// Define functions */
// -------------------------------------------------------------
void ReadFEandPDFiles(double *PDCoords, DENS_VEC &lp_ext_nds_bc2PD,
    DENS_VEC &lp_int_nds_bc2FE,
    DENS_VEC &lp_ext_coo_bc2PD,
    DENS_VEC &lp_int_coo_bc2FE,
    DENS_VEC &lp_ov_nds)
{

  int nodenum, lp_nNds_bc2PD, lp_nNds_bc2FE, lp_nNds_overl;
  int a1, offset;
  double ax, ay;
  char buf[100];

  FILE *File1, *File2, *File3;

  // Open test file
  if((File1 = fopen("dump_t1_GroupOut", "r")) == NULL) {
    cout << "Error in opening dump_t1_GroupOut" << endl;
  }

  // Open test file
  if((File2 = fopen("dump_t1_GroupIn", "r")) == NULL) {
    cout << "Error in opening dump_t1_GroupIn" << endl;
  }

  // Open test file
  if((File3 = fopen("dump_t1_GroupOver", "r")) == NULL) {
    cout << "Error in opening dump_t1_GroupOver" << endl;
  }

  /** Read Number nodes File1 */
  fgets(buf,100,File1);
  fgets(buf,100,File1);
  fscanf(File1, "%i", &a1);
  fgets(buf,100,File1);
  fscanf(File1, "%i", &lp_nNds_bc2PD);
  lp_ext_nds_bc2PD.reset(lp_nNds_bc2PD);
  lp_ext_coo_bc2PD.reset(2*lp_nNds_bc2PD);

  fgets(buf,100,File1);
  fgets(buf,100,File1);
  fgets(buf,100,File1);
  fscanf(File1, "%lf%lf", &ax, &ay);
  fscanf(File1, "%lf%lf", &ax, &ay);
  fscanf(File1, "%lf%lf", &ax, &ay);
  fgets(buf,100,File1);

  for (int inode=0; inode<lp_nNds_bc2PD; inode++)
  {
    fscanf(File1, "%i", &nodenum);
    lp_ext_nds_bc2PD(inode) = nodenum-1;
    offset = 3*(nodenum-1);
    lp_ext_coo_bc2PD(2*inode+0)=PDCoords[offset+0];
    lp_ext_coo_bc2PD(2*inode+1)=PDCoords[offset+1];
  }

  /** Read Number nodes File2 */
  fgets(buf,100,File2);
  fgets(buf,100,File2);
  fscanf(File2, "%i", &a1);
  fgets(buf,100,File2);
  fscanf(File2, "%i", &lp_nNds_bc2FE);
  lp_int_nds_bc2FE.reset(lp_nNds_bc2FE);
  lp_int_coo_bc2FE.reset(2*lp_nNds_bc2FE);

  fgets(buf,100,File2);
  fgets(buf,100,File2);
  fgets(buf,100,File2);
  fscanf(File2, "%lf%lf", &ax, &ay);
  fscanf(File2, "%lf%lf", &ax, &ay);
  fscanf(File2, "%lf%lf", &ax, &ay);
  fgets(buf,100,File2);

  for (int inode=0; inode<lp_nNds_bc2FE; inode++)
  {
    fscanf(File2, "%i", &nodenum);
    lp_int_nds_bc2FE(inode) = nodenum-1;
    offset = 3*(nodenum-1);
    lp_int_coo_bc2FE(2*inode+0)=PDCoords[offset+0];
    lp_int_coo_bc2FE(2*inode+1)=PDCoords[offset+1];
  }

  /** Read Number nodes File3  */
  fgets(buf,100,File3);
  fgets(buf,100,File3);
  fscanf(File3, "%i", &a1);
  fgets(buf,100,File3);
  fscanf(File3, "%i", &lp_nNds_overl);
  lp_ov_nds.reset(lp_nNds_overl);

  fgets(buf,100,File3);
  fgets(buf,100,File3);
  fgets(buf,100,File3);
  fscanf(File3, "%lf%lf", &ax, &ay);
  fscanf(File3, "%lf%lf", &ax, &ay);
  fscanf(File3, "%lf%lf", &ax, &ay);
  fgets(buf,100,File3);

  for (int inode=0; inode<lp_nNds_overl; inode++)
  {
    fscanf(File3, "%i", &nodenum);
    lp_ov_nds(inode) = nodenum-1;
  }

  fclose (File1);
  fclose (File2);
  fclose (File3);
}
