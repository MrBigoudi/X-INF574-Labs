#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

/**
 * A class for dealin with the computation of Bezier curves and their properties
 **/
class Bezier{

	public:
		/**
		* An iterative implementation of the De Casteljau algorithm for computing a Bezier curve
		*
		* @param V  the vertices of the control polygon
		* @param t   an input parameter in [0, 1]
		*
		* @return    the point B(t), obtaining by evaluating the curve for the value 't'
		**/
		MatrixXd de_casteljau(const MatrixXd &V, double t){
			int nV = V.rows();         // number of points in the control polygon
			int degree = V.rows() - 1; // degree of the curve
			
			MatrixXd m = V;
			for(int i=0; i<degree; i++){
				m = de_casteljau_intermediate(m, t);
				//std::cout << "step " << i << " : " << m.rows() << std::endl;
				//std::cout << m << std::endl;
			}
			return m;
		}

		MatrixXd de_casteljau_intermediate(const MatrixXd &V, double t){
			int nV = V.rows();         // number of points in the control polygon
			int degree = V.rows() - 1; // degree of the curve

			MatrixXd m = MatrixXd::Zero(nV-1,3);
			// for each patch
			for(int i=0; i<degree; i++){
				m.row(i) = (1-t)*V.row(i) + t*V.row(i+1);
			}
			return m;
		}

		/**
		* Plot the curve, for t=0, dt, 2*dt, ..., 1, with a given 'resolution' <br>
		* where dt=1/(resolution-1)
		*
		* @param resolution  number of points to be evaluated on the curve
		*/
		MatrixXd plot_curve(const MatrixXd &V, int resolution){
			double dt = 1. / (resolution - 1);
			MatrixXd m = MatrixXd::Zero(resolution,3);
			for(int i=0; i<resolution; i++){
				m.row(i) = plot_curve(V, resolution, i*dt);
			}
			return m;
		}

		MatrixXd plot_curve(const MatrixXd &V, int resolution, double t){
			//double dt = t/ (resolution - 1);
			return de_casteljau(V, t);
		}


		/**
		* Perform the subdivision (once) of the Bezier curve (with parameter t) <br>
		* Return two Bezier curves (with 'n' control points each)
		*/
		vector<MatrixXd> subdivide(const MatrixXd &V, double t){
			//vector<MatrixXd> curves{}; // the result: store the 2 curves obtained after subdivision
			// To be completed (ex 1)
			int dim = V.rows();
			MatrixXd m1 = MatrixXd::Zero(dim, 3);
			MatrixXd m2 = MatrixXd::Zero(dim, 3);

			MatrixXd m = V;
			m1.row(0) = m.row(0);
			m2.row(dim-1) = m.row(dim-1);
			for(int i=0; i<dim-1; i++){
				m = de_casteljau_intermediate(m, t);
				m1.row(i+1) = m.row(0);
				m2.row(dim-1-1-i) = m.row(m.rows()-1);
			}

			return {m1,m2};
		}

		/**
		* Plot the curve using recursive subdivision <br>
		*
		* @param levels  number of levels of subdivisions
		* @return  a polyline representing the curve to be rendered: this is obtained by concantenation of
		* the control points of all subdivided curves
		*/
		MatrixXd subdivision_plot(const MatrixXd &V, int levels){
			std::cout << "computing recursive subdivision " << std::endl;
			// to be completed
			if(levels == 0){
				return {V};
			}
			vector<MatrixXd> sub = subdivide(V, 0.5f);
			// create the matrix from the list
			MatrixXd m1 = subdivision_plot(sub[0], levels-1);
			MatrixXd m2 = subdivision_plot(sub[1], levels-1);

			MatrixXd res(m1.rows()+m2.rows(), m1.cols());
			res << m1, m2;
			return res;
		}

		/**
		* Compute the tangent of a given curve c(t) for a given parameter t0
		*
		* @param V  the vertices of the control polygon
		* @param t0   an input parameter in [0, 1]
		*
		* @return    the tangent at c(t0)
		**/
		MatrixXd compute_tangent(const MatrixXd &V, double t0){
			int n = V.rows(); // number of points in the control polygon
			MatrixXd V_prime = MatrixXd::Zero(n-1, 3);
			for(int i=0; i<n-1; i++){
				V_prime.row(i) = n*(V.row(i+1)-V.row(i));
			}
	
			MatrixXd tan = de_casteljau(V_prime, t0);
			return tan.normalized();
		}

		/**
		 * Compute the cross product of two MatrixXd
		 * @param m1 The first matrix
		 * @param m2 The second matrix
		 * @return The cross product as a MatrixXd
		 * @cond m1 and m2 must be 1x3 matrices
		*/
		MatrixXd getCross(const MatrixXd& m1, const MatrixXd& m2) const {
			MatrixXd res(1,3);
			res << (m1(0,1)*m2(0,2) - m2(0,1)*m1(0,2)),
				   (m1(0,2)*m2(0,0) - m2(0,2)*m1(0,0)),
				   (m1(0,0)*m2(0,1) - m2(0,0)*m1(0,1));
			return res;
		}	

		/**
		* Compute the normal vector of a given curve c(t) for a given parameter t0
		*
		* @param V  the vertices of the control polygon
		* @param t0   an input parameter in [0, 1]
		*
		* @return    the normal at c(t0)
		**/
		MatrixXd compute_normal(const MatrixXd &V, double t0){
			int n = V.rows(); // number of points in the control polygon
			// To be completed
			MatrixXd tan = compute_tangent(V, t0);
			MatrixXd tmp(1,3);
			tmp << 0, 0, 1;
			tmp = getCross(tan, tmp);
			return tmp.normalized();
		}

		/**
		* Compute a loop of points around a curve c(t) for a given parameter t0
		* The points belongs on a circle lying the hyperplane passing through c(t0) and orthogonal to tangent(t0)
		*
		* @param V  the vertices of the control polygon
		* @param t0   an input parameter in [0, 1]
		*
		* @return    a loop of vertices on the hyperplane passing through c(t0) and orthogonal to tangent(t0)
		**/
		MatrixXd compute_loop_of_vertices(const MatrixXd &V, double t0, int k, double radius){
			int n = V.rows(); // number of points in the control polygon
			// To be completed
			MatrixXd tan = compute_tangent(V, t0);
			MatrixXd ct = de_casteljau(V, t0);
			MatrixXd norm = compute_normal(V, t0);
			MatrixXd cross = getCross(tan, norm);

			function<Vector3d(float theta)> circle = [=](float theta){
				float xTmp = ct(0,0) + radius * (cos(theta) * cross(0,0) + sin(theta) * norm(0,0));
				float yTmp = ct(0,1) + radius * (cos(theta) * cross(0,1) + sin(theta) * norm(0,1));
				float zTmp = ct(0,2) + radius * (cos(theta) * cross(0,2) + sin(theta) * norm(0,2));
				return Vector3d(xTmp, yTmp, zTmp);
			};

			const float PI = 3.1416;
			float dt = 2*PI / k;

			MatrixXd res = MatrixXd::Zero(k,3);
			for(int i=0; i<k; i++){
				res.row(i) = circle(i*dt);
			}

			return res;
		}
};
