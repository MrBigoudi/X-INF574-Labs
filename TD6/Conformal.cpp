#include <cmath>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/harmonic.h>
#include "Conformal.h"

/**
 * Given the index of 3 vertices, get the angle between them
 * @param vIdx The vertex origin of the angle
 * @param v1Idx The first vertex
 * @param v2Idx The second vertex
 * @return The angle between the two vertices
*/
float getTipAngle(const MatrixXd &V, int vIdx, int v1Idx, int v2Idx){
	VectorXd v  = V.row(vIdx);
	VectorXd v1 = V.row(v1Idx);
	VectorXd v2 = V.row(v2Idx);

	v1 -= v;
	v2 -= v;

	float dot = v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();

	float res = acos(dot / (v1.norm()*v2.norm()));

	return res;
}

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

	std::vector<Eigen::Triplet<double>> list;

	for(int f=0; f<F.rows(); f++){
		int viIdx = F(f, 0);
		int vjIdx = F(f, 1);
		int vkIdx = F(f, 2);

		float alpha_jk = getTipAngle(V, viIdx, vjIdx, vkIdx);
		float alpha_ki = getTipAngle(V, vjIdx, vkIdx, viIdx);
		float alpha_ij = getTipAngle(V, vkIdx, viIdx, vjIdx);

		float dist_ij = viIdx - vjIdx;
		float dist_jk = vjIdx - vkIdx;
		float dist_ki = vkIdx - viIdx;

		float val_ij = cotangent(alpha_ij)*dist_ij*dist_ij;
		float val_jk = cotangent(alpha_jk)*dist_jk*dist_jk;
		float val_ki = cotangent(alpha_ki)*dist_ki*dist_ki;

		list.push_back({viIdx, vjIdx, val_ij});
		list.push_back({vjIdx, viIdx, val_ij});
		list.push_back({vjIdx, vkIdx, val_jk});
		list.push_back({vkIdx, vjIdx, val_jk});
		list.push_back({vkIdx, viIdx, val_ki});
		list.push_back({viIdx, vkIdx, val_ki});
	}
		
	// convert the matrix to sparse matrix
	Dirichlet.resize(V.rows(), V.rows());
	Dirichlet.setFromTriplets(list.begin(), list.end());
}

void ConformalParametrization::compute_harmonic_parametrization(MatrixXd &V_uv_harmonic){
	igl::boundary_loop(F, b);

	// map boundaries in the unit circle
	MatrixXd bc = MatrixXd::Zero(b.size(), 2);
	for(int i=0; i< b.size(); i++){
		double angle = 2.0 * M_PI * i / b.size();
		bc(i,0) = cos(angle);
		bc(i,1) = sin(angle);
	}

	// Dirichlet*u = B
	VectorXd B = VectorXd::Zero(V.rows(), 1);
	// prepare data for min_quad_with_fixed
	SparseMatrix<double> Q = Dirichlet;

    igl::min_quad_with_fixed_data<double> data;
	printf("test\n\n");
    igl::min_quad_with_fixed_precompute(Q.eval(), b, SparseMatrix<double>(), true, data);
	printf("test\n\n");
    // solve the minimization problem
    VectorXd minimal_z;
    igl::min_quad_with_fixed_solve(data, B, bc, VectorXd(), minimal_z);
	printf("test\n\n");

    // extract the solution and overwrite V_uv_harmonic
    V_uv_harmonic = Map<MatrixXd>(minimal_z.data(), minimal_z.size() / 2, 2);
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



