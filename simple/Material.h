#ifndef MATERIAL_H
#define MATERIAL_H

using namespace std;

/**
 *  @class  Material
 *  @brief  Base class for a
 */
class Material {

public:
	Material();

	Material(double mu, double nu, int ndof, int ncoord);

	double **** material_stiffness(int problem_type) const;

	~Material();

	enum FE_ProblemType { PLANESTRAIN, PLANESTRESS};

protected:

	//
	double mu_;

	//
	double nu_;

	//
	int ndof_;

	//
	int ncoord_;

	//
	double ****C_;

};

#endif // MATERIAL_H
