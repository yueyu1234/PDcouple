#include "Material.h"

#include <iostream>

using namespace std;

// -------------------------------------------------------------
// -------------------------------------------------------------
//   class Material
// -------------------------------------------------------------
// -------------------------------------------------------------

Material::Material()
{

}

Material::Material(double mu_, double nu_, int ndof_, int ncoord_):
mu_(mu_),
nu_(nu_),
ndof_(ndof_),
ncoord_(ncoord_)
{
	// allocate tensor
	C_ = new double***[ndof_];
	for (int i = 0; i < ndof_; ++i) {
		C_[i] = new double**[ncoord_];

		for (int j = 0; j < ncoord_; ++j) {
			C_[i][j] = new double*[ndof_];

			for (int k = 0; k < ndof_; ++k) {
				C_[i][j][k] = new double[ncoord_];
			}
		}
	}
	//C(ndof_,ncoord_,ndof_,ncoord_);
}

Material::~Material()
{
	// de-allocate memory
	for (int i = 0; i < ndof_; ++i) {
		for (int j = 0; j < ncoord_; ++j) {
			for (int k = 0; k < ndof_; ++k)
				delete [] C_[i][j][k];

			delete [] C_[i][j];
		}
		delete [] C_[i];
	}
	delete [] C_;
}

double **** Material::material_stiffness(int problemtype) const
{

//	PLANESTRAIN = 1;
/*	if (ncoord_ == 2) {
		for(int i = 0; i < 2; ++i) {
			for(int j = 0; j < 2; ++j) {
				for(int k = 0; k < 2; ++k) {
					for(int l = 0; l < 2; ++l) {
						if (problemtype == 1) {
							if (i==j && k==l) {
								C_[i][j][k][l] = C_[i][j][k][l] + 2*mu_*nu_/(1-2*nu_);
							}
							else {
								if (i==j && k==l) {
									C_[i][j][k][l] = C_[i][j][k][l] + 2*mu_*nu_/(1-nu_);
								}
							}
							if (i==l && k==j) {
								C_[i][j][k][l] = C_[i][j][k][l] + mu_;
							}
							if (i==k && j==l) {
								C_[i][j][k][l] = C_[i][j][k][l] + mu_;
							}
						}
					}
				}
			}
		}
        }*/
	if (ncoord_ == 2) {
cout << "HERE!!!!!!!!!!!!!!!!!!!!!!ncoords==2" << endl; cout <<endl; cout <<endl;
		for(int i = 0; i < 2; ++i) {
			for(int j = 0; j < 2; ++j) {
				for(int k = 0; k < 2; ++k) {
					for(int l = 0; l < 2; ++l) {
						if (problemtype == 1) {
							if (i==j && k==l) {
                                                                if (i==k)
								    C_[i][j][k][l] = mu_*(1.-nu_)/(1.+nu_)/(1.-2.*nu_);
                                                                else
                                                                    C_[i][j][k][l] = nu_*mu_/(1.+nu_)/(1.-2.*nu_);
							}
                                                        else if (i==l && k==j) {
								C_[i][j][k][l] = mu_/(1.+nu_)/2.;
							}
                                                        else if (i==k && j==l) {
								C_[i][j][k][l] = mu_/(1.+nu_)/2.;
							}
						}
					}
				}
			}
		}
	}
	else if (ncoord_ == 3) {

		for(int i = 0; i < 3; ++i) {
			for(int j = 0; j < 3; ++j) {
				for(int k = 0; k < 3; ++k) {
					for(int l = 0; l < 3; ++l) {
						if (i==j && k==l) {
							C_[i][j][k][l] = C_[i][j][k][l] + 2.*mu_*nu_/(1.-2.*nu_);
						}
						if (i==k && j==l) {
							C_[i][j][k][l] = C_[i][j][k][l] + mu_;
						}
						if (i==l && j==k) {
							C_[i][j][k][l] = C_[i][j][k][l] + mu_;
						}
					}
				}
			}
		}
	}
	return C_;
}
