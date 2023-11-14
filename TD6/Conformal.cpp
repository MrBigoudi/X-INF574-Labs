#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/harmonic.h>
#include "Conformal.h"



ConformalParametrization::ConformalParametrization(const MatrixXd &V0, const MatrixXi &F0){
	V =  V0;
	F = F0;
	igl::boundary_loop(F,boundary);

	//fix points on the boundary -- experiment here by fixing different points
	b.resize(2,1);
	b(0) = boundary(0);
	b(1) = boundary(boundary.size()/2);
	//b(1) = boundary(37);

	//where should we map them

	MatrixXd bc(2,2);
	bc<<0,0,1,0;
}

double ConformalParametrization::cotangent(double x){
	return 1/(tan(x));
}

//compute the matrix for dirichlet energy of the system
void ConformalParametrization::compute_dirichlet(SparseMatrix<double> &Dirichlet){
	
	//TODO
	
}

void ConformalParametrization::compute_harmonic_parametrization(MatrixXd &V_uv_harmonic){
	
	//TODO

}

void ConformalParametrization::compute_area(SparseMatrix<double> &Area){
	//TODO
	
}

void ConformalParametrization::compute_conformal_energy(SparseMatrix<double> &ConformalEnergy){
	//TODO
	
}

void ConformalParametrization::minimize_energy(MatrixXd &V_uv){
	//TODO
}

void ConformalParametrization::build_parametrizations(){
	compute_dirichlet(Dirichlet);
	compute_harmonic_parametrization(V_uv_harmonic);
	compute_area(Area);
	compute_conformal_energy(ConformalEnergy);
	minimize_energy(V_uv);
}

void ConformalParametrization::compute_conformal_factor(MatrixXd &conformal_factor){
	//TODO
}

void ConformalParametrization::compute_conformal_coloring(const MatrixXd &conformal_factor, MatrixXd Color){
	//TODO
}



