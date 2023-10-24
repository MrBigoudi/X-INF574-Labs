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

#include "HalfedgeBuilder.cpp"


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
		int nEdge = he.getOpposite(he.getNext(beg));
		while(nEdge != beg){
			nEdge = he.getOpposite(he.getNext(nEdge));
			res++;
		}
		return res;
	}

	/**
	 * Return the degree of a given vertex 'v'
	 * Turn around vertex 'v' in CCW order
	 **/
	int vertexDegreeCCW(HalfedgeDS he, int v)
	{
	        // TO BE COMPLETED
	        
	        return 0;
	}

	/**
	 * Compute the vertex normals (he)
	 * Exercice
	 **/
	void vertexNormals(HalfedgeDS he) {
		// TODO --fill the matrix he_N_vertices
		
	}

	// Compute lib_per-vertex normals
	/**
	 * Compute the vertex normals (global, using libiGl data structure)
	 * Exercice
	**/
	void lib_vertexNormals() {
		// TODO -- fill the matrix lib_N_vertices
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
		lib_Deg_vertices = MatrixXi(V.rows(),1);
		for(int i=0; i<V.rows(); i++){
			lib_Deg_vertices(i,0) = V.row(i).size();
		}
		auto end = std::chrono::high_resolution_clock::now(); // for measuring time performances
		std::cout << "Time spent: " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << "ms" << std::endl;
	}
    	
	/**
	* Return the number of boundaries of the mesh
	* Exercice
	**/
	int countBoundaries(HalfedgeDS he){
		// TO BE COMPLETED
		return 0;
	}

	/**
	 * @brief Computes the dual voronoi areas of every vertex. The MatrixXd voronoi has to be filled.
	 * 
	 * @return ** void 
	 */
	void compute_voronoi_area(){
		//TODO
	}

	void compute_gaussian_curvature(HalfedgeDS he){
		// TODO

	}

	void compute_lib_gaussian_curvature(){
		// TODO
	}

	void updateColorGaussianCurvature(MatrixXd &C, double t = 0.0){
		//TODO 
	}

}


namespace laplacian{

	double cotangent(double x){
		return 1/(tan(x));
	}

	void buildLaplacian(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double> &L, SparseMatrix<double> &A, SparseMatrix<double> &Ainv, SparseMatrix<double> &Delta){

		//fill the Dirichlet Matrix

		//build the matrices L, A, A_inv, Delta as seen in the Lecture 
		MatrixXd L_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd A_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd A_inv_dense = MatrixXd::Zero(V.rows(),V.rows());
		MatrixXd Delta_dense = MatrixXd::Zero(V.rows(),V.rows());

		// TODO -> Fill the matrices
		
		
		// Now convert the matrices to sparse matrices

		L = L_dense.sparseView();
		A = A_dense.sparseView();
		Ainv = A_inv_dense.sparseView();
		Delta = Delta_dense.sparseView();
	}

	void setHeat(MatrixXd & u){
		//TODO
	}

	void heat_step_explicit(const SparseMatrix<double> &Delta, MatrixXd & u, double time_step ){
		//TODO
	}

	void heat_step_implicit(const SparseMatrix<double> &L,const SparseMatrix<double> &A, MatrixXd & u, double time_step){
		//TODO 
	}

	void updateColorHeat(const MatrixXd & u, MatrixXd &C){
		//TODO
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

    igl::readOFF("./data/cat0.off",V,F); // change this line depending on your system

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

	// TASK 1.3
	// Compute normals

	// Use these pre defined methods as a check for your implementation 
	// Compute per-face normals
	igl::per_face_normals(V,F,N_faces);
	// Compute per-vertex normals
	igl::per_vertex_normals(V,F,N_vertices);
	
	// Compute lib_per-vertex normals
	// mesh_methods::lib_vertexNormals();

	// Compute he_per-vertex normals
    // mesh_methods::vertexNormals(he); 

	// TASK 1.4 
	// compute number of boundaries
	// int B=mesh_methods::countBoundaries(he);  
	// std::cout << "The mesh has " << B << " boundaries" << std::endl;




	auto time = (std::chrono::high_resolution_clock::now());
	std::cout<<
    "Press '1' for per-face normals calling pre-defined functions of LibiGL."<<std::endl<<
	"Press '2' for per-vertex normals calling pre-defined functions of LibiGL."<<std::endl<<
    "Press '3' for lib_per-vertex normals using face-vertex structure of LibiGL ."<<std::endl<<
	"Press '4' for HE_per-vertex normals using HalfEdge structure."<<std::endl;

	
	// TASK 2
	// compute the gaussian curvature
	// mesh_methods::compute_voronoi_area();
	// mesh_methods::compute_gaussian_curvature(he); //here use the half edge structure
	// mesh_methods::compute_lib_gaussian_curvature(); //here use the libigl face based datastructure
	
	
	
	
	
	

	/*---------------- Uncomment this part for the build of the Laplacian ----------------*/
	// TASK 3
	// laplacian::buildLaplacian(V,F,L,A,Ainv,Delta);
	// MatrixXd u = MatrixXd::Zero(V.rows(),1);
	// MatrixXd u_prev = MatrixXd::Zero(V.rows(),1);
	// laplacian::setHeat(u);
	// laplacian::updateColorHeat(u, C);
	// viewer.data().set_colors(C);


	/*---------------- Uncomment this part for the animation ----------------*/

	
	viewer.core().is_animating = true;  // Enable animation
    double t = 0.0;  // Initialize time
	double time_step = 0.1; // experiment with this parameter

	// std::this_thread::sleep_for(std::chrono::milliseconds(10000));
    viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool {
	    // Use a delay to control the speed of the color change
	    // std::this_thread::sleep_for(std::chrono::milliseconds(5));
	    // Calculate new vertex colors based on your time-dependent function
		
		//--------- for time dependent coloring ---------
		// TASK 1.2
		// mesh_methods::updateColor(C,t);
		// viewer.data().set_colors(C);
		
		//--------- for the heat flow ---------
		// laplacian::setHeat(u);	
		// laplacian::heat_step_explicit(Delta,u,time_step);
		// laplacian::heat_step_implicit(L,A,u,time_step);

		// laplacian::updateColorHeat(u,C);
		// viewer.data(0).set_colors(C);
		// viewer.data().set_colors(C);

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
