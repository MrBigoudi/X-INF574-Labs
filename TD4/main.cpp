#include <cmath>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <vector>

#include "HalfedgeBuilder.cpp"
#include "igl/PI.h"
#include "igl/cross.h"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;


MatrixXd N_faces;   //computed calling pre-defined functions of LibiGL
MatrixXd N_vertices; //computed calling pre-defined functions of LibiGL


MatrixXd lib_N_vertices;  //computed using face-vertex structure of LibiGL
MatrixXi lib_Deg_vertices;//computed using face-vertex structure of LibiGL

MatrixXd he_N_vertices; //computed using the HalfEdge data structure

MatrixXd curvatures; //computed using the HalfEdge data structure
MatrixXd lib_curvatures; // computed using the face-vertex structure of LibiGL

MatrixXd voronoi; // area of the dual voronoi areas of every vertex


// colors for the dfs
enum COLORS{
	WHITE,
	BLACK,
	GREY,
};

SparseMatrix<double> L; //cotangent
SparseMatrix<double> A; //normalization for the laplacian
SparseMatrix<double> Ainv; //inverse of A for performance
SparseMatrix<double> Delta; //laplacian

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
	switch(key){
		case '1':
			viewer.data().set_normals(N_faces);
			return true;
		case '2':
			viewer.data().set_normals(N_vertices);
			return true;
		case '3':
			viewer.data().set_normals(lib_N_vertices);
			return true;
		case '4':
			viewer.data().set_normals(he_N_vertices);
			return true;
		default: break;
	}
	return false;
}

double cotangent(double x){
	return 1/(tan(x));
}

namespace mesh_methods{
	
	void updateColor(MatrixXd &C, double t = 0.0){
		// TODO -- TASK 1.2
		C *= cos(t) + sin(t);
	}

	/**
	 * Return the degree of a given vertex 'v'
	 * Exercice
	 **/
	int vertexDegree(HalfedgeDS he, int v) {
		// TO BE COMPLETED
		int res = 0;
		int beg = he.getEdge(v);
		int next = beg;
		do{
			next = he.getOpposite(he.getNext(next));
			res++;
		}while(next != beg);
		return res;
	}

	/**
	 * Return the degree of a given vertex 'v'
	 * Turn around vertex 'v' in CCW order
	 **/
	int vertexDegreeCCW(HalfedgeDS he, int v){
		// TO BE COMPLETED
		return 0;
	}

	/**
	 * Compute the vertex normals (he)
	 * Exercice
	 **/
	void vertexNormals(HalfedgeDS he) {
		// TODO --fill the matrix he_N_vertices
		int nbVert = he.sizeOfVertices();
		he_N_vertices = MatrixXd::Zero(nbVert, 3);
		for(int v=0; v<nbVert; v++){
			MatrixXd normal = MatrixXd::Zero(1,3);
			int beg = he.getEdge(v);
			int prev = beg;
			int next = he.getOpposite(he.getNext(prev));
			// find all the faces surrounding it
			do{
				// get the face normal
				MatrixXd e1 = MatrixXd::Zero(1,3);
				MatrixXd e2 = MatrixXd::Zero(1,3);
				int v1 = he.getTarget(he.getOpposite(prev));
				int v2 = he.getTarget(he.getOpposite(next));

				int tmp = v;
				tmp = min(tmp, v1);


				e1.row(0) = V.row(v1) - V.row(v);
				e2.row(0) = V.row(v2) - V.row(v);
				MatrixXd normalTmp = MatrixXd::Zero(1,3);
				igl::cross(e2, e1, normalTmp);
				normalTmp.normalize();
				normal += normalTmp;

				prev = next;
				next = he.getOpposite(he.getNext(next));
			}while(prev != beg);
			he_N_vertices.row(v) = normal.normalized();
		}
	}

	// Compute lib_per-vertex normals
	/**
	 * Compute the vertex normals (global, using libiGl data structure)
	 * Exercice
	**/
	void lib_vertexNormals() {
		// TODO -- fill the matrix lib_N_vertices
		lib_N_vertices = MatrixXd::Zero(V.rows(), 3);
		// get the faces normals
		for(int f=0; f<F.rows(); f++){
			// get the vertices
			int v1 = F.row(f).x();
			int v2 = F.row(f).y();
			int v3 = F.row(f).z();

			// get the face normal
			MatrixXd e1 = MatrixXd::Zero(1,3);
			MatrixXd e2 = MatrixXd::Zero(1,3);
			e1.row(0) = V.row(v3) - V.row(v2);
			e2.row(0) = V.row(v1) - V.row(v2);
			MatrixXd normal = MatrixXd::Zero(1,3);
			igl::cross(e1, e2, normal);
			normal.normalize();

			// add the normal to the vertices normals
			lib_N_vertices.row(v1) += normal;
			lib_N_vertices.row(v2) += normal;
			lib_N_vertices.row(v3) += normal;
		}

		// normalize the vertices normals
		for(int v=0; v<V.rows(); v++){
			lib_N_vertices.row(v).normalize();
		}
	}

	/**
	 * Return the number of occurrence of vertex degrees: for d=3..n-1
   	 * Exercice
	 **/
	void vertexDegreeStatistics(HalfedgeDS he) {
		std::cout << "Computing vertex degree distribution...";
		auto start = std::chrono::high_resolution_clock::now(); // for measuring time performances
		int n = he.sizeOfVertices();
		std::map<int, int> nbDegrees;
		for(int v=0; v<n; v++){
			int deg = vertexDegree(he, v);
			if(nbDegrees.find(deg) != nbDegrees.end())
				nbDegrees[deg]++;
			else
				nbDegrees[deg] = 1;
		}
		for(int d=3; d<n; d++){
			int deg = 0;
			if(nbDegrees.find(d-3) != nbDegrees.end())
				deg = nbDegrees[d-3];
	    	//std::cout << "number of degrees " << d << " = "<< deg <<std::endl;
		}
		auto end = std::chrono::high_resolution_clock::now(); // for measuring time performances
		std::cout << "Time spent: " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << "ms" << std::endl;
	}

	// Compute lib_vertex degrees
	/**
	 *  (global, using libiGl data structure)
	 * 
	 * 
	**/
	void lib_vertexDegrees() {
		//TODO -- fill the matrix lib_Deg_vertices
		std::cout << "Computing vertex degrees using libigl...";
		auto start = std::chrono::high_resolution_clock::now(); // for measuring time performances
		lib_Deg_vertices = MatrixXi::Zero(V.rows(),1);

		for(int f=0; f<F.rows(); f++){
			for(int v=0; v<F.row(f).size(); v++){
				lib_Deg_vertices(F.row(f)[v], 0)++;
				// lib_Deg_vertices(v1,0)+=2;
				// lib_Deg_vertices(v2,0)+=2;
				// lib_Deg_vertices(v3,0)+=2;
			}
		}
		// do not need to divide by 2 if only adding 1
		// for(int v=0; v<V.rows(); v++){
		// 	lib_Deg_vertices(v,0) /= 2;
		// }
		auto end = std::chrono::high_resolution_clock::now(); // for measuring time performances
		std::cout << "Time spent: " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << "ms" << std::endl;
	}

	/**
	 * Get all the neighbours of a given vertex
	 * @param he The half edge structure
	 * @param v The current vertex
	 * @return All its neighbours as a vector
	*/
	std::vector<int> getAllVerticesNeighbours(HalfedgeDS& he, int v){
		std::vector<int> neighbours = std::vector<int>();
		int beg = he.getEdge(v);
		int next = beg;
		do{
			neighbours.push_back(he.getTarget(he.getOpposite(next)));
			next = he.getOpposite(he.getNext(next));
		}while(next != beg);

		return neighbours;
	}

	/**
	 * Do an iteration of the DFS visit
	 * @param he The half edge structure
	 * @param e The current edge
	 * @param visited The state of all the vertices
	*/
	void dfsVisit(HalfedgeDS& he, int e, MatrixXi& visited){
		visited.row(e)[0] = GREY;
		int neighbour = -1;
		neighbour = he.getNext(e);
		if(visited.row(neighbour)[0] == WHITE){
			dfsVisit(he, neighbour, visited);
		}
		visited.row(e)[0] = BLACK;
	}

	/**
	 * Tells if an edge is at a boundary
	 * @param he The half edge structure
	 * @param e The edge to check
	*/
	bool isAtBoundary(HalfedgeDS he, int e){
		return (e != he.getNext(he.getNext(he.getNext(e))));
	}
	
    	
	/**
	* Return the number of boundaries of the mesh
	* Exercice
	**/
	int countBoundaries(HalfedgeDS he){
		// TO BE COMPLETED
		std::cout << "he = " << he.sizeOfHalfedges() << ", nbFaces = " << F.rows() << std::endl;
		if(he.sizeOfHalfedges() == 3*F.rows()) return 0;

		int nbConnectedComponents = 0;
		MatrixXi visited = MatrixXi::Zero(he.sizeOfHalfedges(),1);
		for(int e=0; e<visited.rows(); e++){
			visited.row(e)[0] = WHITE;
		}

		// do a DFS
		for(int e=0; e<visited.rows(); e++){
			if(visited.row(e)[0] == WHITE){
				// left of the boundary
				if(isAtBoundary(he, e) && isAtBoundary(he, he.getNext(e))){
					dfsVisit(he, e, visited);
					nbConnectedComponents++;
					continue;
				}
				// right of the boundary
				if(isAtBoundary(he, e) && isAtBoundary(he, he.getPrev(he.getOpposite(e)))){
					visited.row(e)[0] = BLACK;
					dfsVisit(he, he.getOpposite(e), visited);
					nbConnectedComponents++;
					continue;
				}
			}
		}

		return nbConnectedComponents;
	}

	void updateColorBoundaries(MatrixXd &C, HalfedgeDS& he){
		//TODO 
		for(int e=0; e<he.sizeOfHalfedges(); e++){
			MatrixXd color = MatrixXd(1, 3);
			if(isAtBoundary(he, e)){
				color << 1,0,0;
			} else {
				color << 1, 1, 1;
			}
			C.row(he.getTarget(e)) = color;
			C.row(he.getTarget(he.getOpposite(e))) = color;
		}
	}

	/**
	 * Given the index of 3 vertices, get the angle between them
	 * @param vIdx The vertex origin of the angle
	 * @param v1Idx The first vertex
	 * @param v2Idx The second vertex
	 * @return The angle between the two vertices
	*/
	float getTipAngle(int vIdx, int v1Idx, int v2Idx){
		VectorXd v  = V.row(vIdx);
		VectorXd v1 = V.row(v1Idx);
		VectorXd v2 = V.row(v2Idx);

		v1 -= v;
		v2 -= v;

		float dot = v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();

		float res = acos(dot / (v1.norm()*v2.norm()));

		return res;
	}

	/**
	 * @brief Computes the dual voronoi areas of every vertex. The MatrixXd voronoi has to be filled.
	 * 
	 * @return ** void 
	 */
	void compute_voronoi_area(){
		//TODO
		voronoi = MatrixXd::Zero(V.rows(),1);
		for(int f=0; f<F.rows(); f++){
			int v1Idx = F.row(f)[0];
			int v2Idx = F.row(f)[1];
			int v3Idx = F.row(f)[2];

			VectorXd v1 = V.row(v1Idx);
			VectorXd v2 = V.row(v2Idx);
			VectorXd v3 = V.row(v3Idx);

			float dist_12 = (v1-v2).norm()*(v1-v2).norm();
			float dist_23 = (v2-v3).norm()*(v2-v3).norm();
			float dist_31 = (v3-v1).norm()*(v3-v1).norm();

			float alpha_1 = getTipAngle(v1Idx, v2Idx, v3Idx);
			float alpha_2 = getTipAngle(v2Idx, v3Idx, v1Idx);
			float alpha_3 = getTipAngle(v3Idx, v1Idx, v2Idx);

			if(alpha_1 < igl::PI / 2.0f){
				voronoi.row(v2Idx)[0] += cotangent(alpha_1)*dist_23;
				voronoi.row(v3Idx)[0] += cotangent(alpha_1)*dist_23;
			}

			if(alpha_2 < igl::PI / 2.0f){
				voronoi.row(v1Idx)[0] += cotangent(alpha_2)*dist_31;
				voronoi.row(v3Idx)[0] += cotangent(alpha_2)*dist_31;
			}

			if(alpha_3 < igl::PI / 2.0f){
				voronoi.row(v1Idx)[0] += cotangent(alpha_3)*dist_12;
				voronoi.row(v2Idx)[0] += cotangent(alpha_3)*dist_12;
			}
		}

		for(int v=0; v<V.rows(); v++){
			voronoi.row(v)[0] /= 8.;
			//std::cout << "Voronoi[" << v << "] = " << voronoi.row(v)[0] << std::endl;
		}
	}

	void compute_gaussian_curvature(HalfedgeDS he){
		// TODO
		curvatures = MatrixXd::Zero(V.rows(),1);
		// for each vertex, get the sumof tip angles
		for(int v=0; v<V.rows(); v++){
			float sum = 0.0f;
			std::vector<int> neighbours = getAllVerticesNeighbours(he, v);
			for(int i=0; i<neighbours.size(); i++){
				int prev = neighbours[i];
				int next = neighbours[(i+1)%neighbours.size()];
				sum += getTipAngle(v, prev, next);
			}
			// get the gaussian curvature
			curvatures.row(v)[0] = (2*igl::PI - sum) / voronoi.row(v)[0];
			//curvatures.row(v)[0] = (2*igl::PI - sum); // libigl doesn't divide by voronoi areas
		}
	}

	void compute_lib_gaussian_curvature(){
		// TODO
		lib_curvatures = MatrixXd::Zero(V.rows(),1);
		// for each face, get the sumof tip angles
		for(int f=0; f<F.rows(); f++){
			int v1Idx = F.row(f)[0];
			int v2Idx = F.row(f)[1];
			int v3Idx = F.row(f)[2];

			VectorXd v1 = V.row(v1Idx);
			VectorXd v2 = V.row(v2Idx);
			VectorXd v3 = V.row(v3Idx);

			lib_curvatures.row(v1Idx)[0] += getTipAngle(v1Idx, v2Idx, v3Idx);
			lib_curvatures.row(v2Idx)[0] += getTipAngle(v2Idx, v3Idx, v1Idx);
			lib_curvatures.row(v3Idx)[0] += getTipAngle(v3Idx, v1Idx, v2Idx);
		}
		// rearrange the curvature for every vertex
		for(int v=0; v<V.rows(); v++){
			lib_curvatures.row(v)[0] = (2*igl::PI - lib_curvatures.row(v)[0]) / voronoi.row(v)[0];
			//lib_curvatures.row(v)[0] = (2*igl::PI - lib_curvatures.row(v)[0]); // libigl doesn't divide by voronoi areas
		}
	}

	void updateColorGaussianCurvature(MatrixXd &C, double t = 0.0){
		//TODO 
		for(int v=0; v<C.rows(); v++){
			MatrixXd color = MatrixXd(1, 3);
			float c = (lib_curvatures.row(v)[0] + 2*igl::PI) / (4*igl::PI);
			// std::cout << "c: " << c << std::endl;
			color << c,0,1-c;
			C.row(v) = color;
		}
	}

}


namespace laplacian{

	void printSparseMatrix(const SparseMatrix<double> &mat, const std::string& name) {
		for (int k = 0; k < mat.outerSize(); ++k) {
			for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
				std::cout << name << "(" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;
			}
		}
	}

	void printMatrixXd(const Eigen::MatrixXd &mat, const std::string& name) {
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				std::cout << name << "(" << i << "," << j << ") = " << mat(i, j) << std::endl;
			}
		}
	}

	void buildLaplacian(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double> &L, SparseMatrix<double> &A, SparseMatrix<double> &Ainv, SparseMatrix<double> &Delta){

		//fill the Dirichlet Matrix

		//build the matrices L, A, A_inv, Delta as seen in the Lecture 
		MatrixXd L_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd A_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd A_inv_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd Delta_dense = MatrixXd::Zero(V.rows(),V.rows());

		// TODO -> Fill the matrices
		for(int i=0; i<V.rows(); i++){
			A_dense(i,i) = voronoi.row(i)[0];
			A_inv_dense(i,i) = 1.0/A_dense(i,i);
		}

		MatrixXd sum = MatrixXd::Zero(V.rows(), 1);
		for(int f=0; f<F.rows(); f++){
			int v1Idx = F.row(f)[0];
			int v2Idx = F.row(f)[1];
			int v3Idx = F.row(f)[2];

			float alpha_1 = mesh_methods::getTipAngle(v1Idx, v2Idx, v3Idx);
			float alpha_2 = mesh_methods::getTipAngle(v2Idx, v3Idx, v1Idx);
			float alpha_3 = mesh_methods::getTipAngle(v3Idx, v1Idx, v2Idx);

			float val3 = (1.0/2.0)*cotangent(alpha_3);
			float val2 = (1.0/2.0)*cotangent(alpha_2);
			float val1 = (1.0/2.0)*cotangent(alpha_1);

			L_dense(v1Idx, v2Idx) += val3;
			Delta_dense(v1Idx, v2Idx) += A_inv_dense(v1Idx,v1Idx) * val3;
			L_dense(v2Idx, v1Idx) += val3;
			Delta_dense(v2Idx, v1Idx) += A_inv_dense(v2Idx,v2Idx) * val3;
			L_dense(v1Idx, v3Idx) += val2;
			Delta_dense(v1Idx, v3Idx) += A_inv_dense(v1Idx,v1Idx) * val2;
			L_dense(v3Idx, v1Idx) += val2;
			Delta_dense(v3Idx, v1Idx) += A_inv_dense(v3Idx,v3Idx) * val2;
			L_dense(v3Idx, v2Idx) += val1;
			Delta_dense(v3Idx, v2Idx) += A_inv_dense(v3Idx,v3Idx) * val1;
			L_dense(v2Idx, v3Idx) += val1;
			Delta_dense(v2Idx, v3Idx) += A_inv_dense(v2Idx,v2Idx) * val1;

			float lDense1 = - val3 - val2;
			float lDense2 = - val3 - val1;
			float lDense3 = - val1 - val2;
			
			L_dense(v1Idx,v1Idx) += lDense1;
			L_dense(v2Idx,v2Idx) += lDense2;
			L_dense(v3Idx,v3Idx) += lDense3;
			
			Delta_dense(v1Idx,v1Idx) += A_inv_dense(v1Idx,v1Idx)*lDense1;
			Delta_dense(v2Idx,v2Idx) += A_inv_dense(v2Idx,v2Idx)*lDense2;
			Delta_dense(v3Idx,v3Idx) += A_inv_dense(v3Idx,v3Idx)*lDense3;
		}
		
		// Now convert the matrices to sparse matrices

		L = L_dense.sparseView();
		A = A_dense.sparseView();
		Ainv = A_inv_dense.sparseView();
		Delta = Delta_dense.sparseView();

		// test
		// MatrixXd HN = MatrixXd::Zero(V.rows(), V.rows());
		// SparseMatrix<double> L2,M,Minv;
		// igl::cotmatrix(V,F,L2);
		// igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
		// igl::invert_diag(M,Minv);
		// HN = Minv*L2;

		// printSparseMatrix(L, "L");
		// printSparseMatrix(L2, "L2");
		// std::cout << std::endl;
		// printSparseMatrix(A, "A");
		// printSparseMatrix(M, "M");
		// std::cout << std::endl;
		// printSparseMatrix(Ainv, "Ainv");
		// printSparseMatrix(Minv, "Minv");
		// std::cout << std::endl;
		// printSparseMatrix(Delta, "Delta");
		// printMatrixXd(HN, "HN");
	}

	void setHeat(MatrixXd & u){
		//TODO
		float minX = 0.0f;
		float minY = 0.0f;
		float maxX = 0.0f;
		float maxY = 0.0f;
		float minZ = 0.0f;
		float maxZ = 0.0f;
		int argminX = 0;
		int argminY = 0;
		int argminZ = 0;
		int argmaxX = 0;
		int argmaxY = 0;
		int argmaxZ = 0;

		float heatValue = 1000.0f;

		for(int v=0; v<V.rows(); v++){
			float valX = V.row(v).x();
			float valY = V.row(v).y();
			float valZ = V.row(v).z();

			if(valX < minX){
				minX = valX;
				argminX = v;
			}
			if(valX > maxX){
				maxX = valX;
				argmaxX = v;
			}
			if(valY < minY){
				minY = valY;
				argminY = v;
			}
			if(valY > maxY){
				maxY = valY;
				argmaxY = v;
			}
			if(valZ < minZ){
				minZ = valZ;
				argminZ = v;
			}
			if(valZ > maxZ){
				maxZ = valZ;
				argmaxZ = v;
			}
		}
		u.row(argminX)[0] = heatValue;
		u.row(argminY)[0] = heatValue;
		u.row(argminZ)[0] = heatValue;
		u.row(argmaxX)[0] = heatValue;
		u.row(argmaxY)[0] = heatValue;
		u.row(argmaxZ)[0] = heatValue;
	}

	void heat_step_explicit(const SparseMatrix<double> &Delta, MatrixXd & u, double time_step ){
		//TODO
		u += time_step*Delta*u;
	}

	void heat_step_implicit(const SparseMatrix<double> &L,const SparseMatrix<double> &A, MatrixXd & u, double time_step){
		//TODO 
		MatrixXd tmp = (A-time_step*L);
		tmp.inverse();
		u = tmp*A*u;
	}

	void updateColorHeat(const MatrixXd & u, MatrixXd &C){
		//TODO
		float heatValue = 1000.0f;
		for(int v=0; v<V.rows(); v++){
			MatrixXd color = MatrixXd(1, 3);
			float c = u.row(v)[0];
			color << c/heatValue,0,0;
			C.row(v) = color;
		}
	}

	void initializeWave(const MatrixXd &V, MatrixXd &u, MatrixXd &u_prev){
		// TODO
	}

	void wave_step(const SparseMatrix<double> &Delta, MatrixXd & u, MatrixXd & u_prev, double time_step ){
		//TODO
	}

	void updateColorWave(const MatrixXd & u, MatrixXd &C){
		//TODO
	}

	void updateMeshPositions(const MatrixXd & u, MatrixXd &V){
		//TODO
	}
}






	


// ------------ main program ----------------
int main(int argc, char *argv[]) {

    // igl::readOFF("./data/bunny_coarser_2.off",V,F); // change this line depending on your system
    igl::readOFF("./data/bunny_coarser.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/bunny_cut.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/bunny_fine.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/bunny_new.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/bunny.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/cat0.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/chandelier.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/cube_open.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/cube_tri.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/face.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/high_genus.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/homer.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/nefertiti.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/sphere.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/star.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/torus.off",V,F); // change this line depending on your system
    // igl::readOFF("./data/venus.off",V,F); // change this line depending on your system

	//print the number of mesh elements
    std::cout << "Points: " << V.rows() << std::endl;


    HalfedgeBuilder* builder=new HalfedgeBuilder();  //

    HalfedgeDS he=builder->createMesh(V.rows(), F);  //


	// Plot the mesh with pseudocolors
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	viewer.callback_key_down = &key_down;
    viewer.data(0).set_mesh(V, F);


	// TASK 1.1

	VectorXd Z;
	Z.setZero(V.rows(),1);
	// Z colors
	// Use the z coordinate as a scalar field over the surface
	// Z= ?
	for(int i=0; i<V.rows(); i++){
		Z[i] = V.row(i).z();
	} 
	MatrixXd C;
	igl::jet(Z,true,C);
	viewer.data().set_colors(C);

	// compute vertex degrees
	mesh_methods::vertexDegreeStatistics(he); 
	mesh_methods::lib_vertexDegrees();
	// for(int v=0; v<V.rows(); v++){
	// 	int deg1 = mesh_methods::vertexDegree(he, v);
	// 	int deg2 = lib_Deg_vertices.row(v).x();
	// 	if(deg1 != deg2){
	// 		std::cout << "wrong output: " << deg1 << " != " << deg2 << std::endl;
	// 	}
	// }

	// TASK 1.3
	// Compute normals

	// Use these pre defined methods as a check for your implementation 
	// Compute per-face normals
	igl::per_face_normals(V,F,N_faces);
	// Compute per-vertex normals
	igl::per_vertex_normals(V,F,N_vertices);
	
	// Compute lib_per-vertex normals
	mesh_methods::lib_vertexNormals();

	// Compute he_per-vertex normals
    mesh_methods::vertexNormals(he); 

	// for(int v=0; v<V.rows(); v++){
	// 	std::cout << "N_vertices[" << v << "] = " << N_vertices.row(v).x() << ", " << N_vertices.row(v).y() << ", " << N_vertices.row(v).z() << std::endl; 
	// 	std::cout << "lib_N_vertices[" << v << "] = " << lib_N_vertices.row(v).x() << ", " << lib_N_vertices.row(v).y() << ", " << lib_N_vertices.row(v).z() << std::endl; 
	// 	std::cout << "he_N_vertices[" << v << "] = " << he_N_vertices.row(v).x() << ", " << he_N_vertices.row(v).y() << ", " << he_N_vertices.row(v).z() << std::endl; 
	// 	std::cout << std::endl;
	// }

	// TASK 1.4 
	// compute number of boundaries
	int B=mesh_methods::countBoundaries(he);  
	std::cout << "The mesh has " << B << " boundaries" << std::endl;
	// mesh_methods::updateColorBoundaries(C, he);
	// viewer.data().set_colors(C);




	auto time = (std::chrono::high_resolution_clock::now());
	std::cout<<
    "Press '1' for per-face normals calling pre-defined functions of LibiGL."<<std::endl<<
	"Press '2' for per-vertex normals calling pre-defined functions of LibiGL."<<std::endl<<
    "Press '3' for lib_per-vertex normals using face-vertex structure of LibiGL ."<<std::endl<<
	"Press '4' for HE_per-vertex normals using HalfEdge structure."<<std::endl;

	
	// TASK 2
	// compute the gaussian curvature
	mesh_methods::compute_voronoi_area();
	mesh_methods::compute_gaussian_curvature(he); //here use the half edge structure
	mesh_methods::compute_lib_gaussian_curvature(); //here use the libigl face based datastructure
	MatrixXd K = MatrixXd::Zero(V.rows(), 1);
	igl::gaussian_curvature(V, F, K);
	// for(int v=0; v<V.rows(); v++){
	// 	std::cout << "he_curv["  << v << "] = " << curvatures.row(v).x() << std::endl; 
	// 	std::cout << "lib_curv[" << v << "] = " << lib_curvatures.row(v).x() << std::endl; 
	// 	std::cout << "igl_curv[" << v << "] = " << K.row(v).x() << std::endl; 
	// 	std::cout << std::endl;
	// }
	// mesh_methods::updateColorGaussianCurvature(C,0.0f);
	// viewer.data().set_colors(C);
	

	/*---------------- Uncomment this part for the build of the Laplacian ----------------*/
	// TASK 3
	laplacian::buildLaplacian(V,F,L,A,Ainv,Delta);
	MatrixXd u = MatrixXd::Zero(V.rows(),1);
	MatrixXd u_prev = MatrixXd::Zero(V.rows(),1);
	laplacian::setHeat(u);
	laplacian::updateColorHeat(u, C);
	viewer.data().set_colors(C);


	/*---------------- Uncomment this part for the animation ----------------*/

	
	viewer.core().is_animating = true;  // Enable animation
    double t = 0.0;  // Initialize time
	// double time_step = 0.1; // experiment with this parameter
	double time_step = 0.00005;
	// double time_step = 1.0;

	// std::this_thread::sleep_for(std::chrono::milliseconds(10000));
    viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool {
	    // Use a delay to control the speed of the color change
	    std::this_thread::sleep_for(std::chrono::milliseconds(5));
	    // Calculate new vertex colors based on your time-dependent function
		
		//--------- for time dependent coloring ---------
		// TASK 1.2
		// mesh_methods::updateColor(C,t);
		// viewer.data().set_colors(C);
		
		//--------- for the heat flow ---------
		laplacian::setHeat(u);	
		laplacian::heat_step_explicit(Delta,u,time_step);
		// laplacian::heat_step_implicit(L,A,u,time_step);

		laplacian::updateColorHeat(u,C);
		//viewer.data(0).set_colors(C);
		viewer.data().set_colors(C);

		//--------- for the wave equation ---------
		// laplacian::wave_step(Delta,u,u_prev,time_step);
		// laplacian::updateColorWave(u,C);
		// laplacian::updateMeshPositions(u,V);
		// viewer.data().clear();
		// viewer.data().set_mesh(V,F);
		// viewer.data().set_colors(C);

		
        
		t += time_step;
        if (t > 1.0) {
            t = 0.0;  // Reset time to loop the animation
        }
	
	    return false;
    };
	

	
	viewer.launch(); // run the viewer
}
