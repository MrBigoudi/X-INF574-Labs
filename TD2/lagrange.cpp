#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include <cmath>
using namespace Eigen;

class LagrangeInterpolation{
  MatrixXd W;
  MatrixXd a;

  MatrixXd Vandermonde(const MatrixXd &V, MatrixXd &W){
    // complete here
    MatrixXd newW = W;
    int len = W.rows();
    for(int i=0; i<len; i++){
      for(int j=0; j<len; j++){
        newW(i,j) = std::pow(V.row(i).x(), j);
      }
    }
    return newW;
  }

public:
  LagrangeInterpolation(const MatrixXd &V1){
    W = MatrixXd::Ones(V1.rows(),V1.rows());
    W = Vandermonde(V1, W);
    // complete here
    //std::cout << W << "\n" << std::endl;
    // create the vector of y
    MatrixXd y = MatrixXd::Ones(V1.rows(), 1);
    for(int i=0; i<V1.rows(); i++){
      y(i,0) = V1.row(i).y();
    }
    //std::cout << y << "\n" << std::endl;
    // create a
    a = W.inverse()*y;
    //std::cout << a << "\n" << std::endl;
  }

  // evaluate interpolated function at time t
  float eval_function(float t){
    // complete here
    // create the phi
    MatrixXd phi = MatrixXd::Ones(W.rows(), 1);
    //std::cout << phi << "\n" << std::endl;
    for(int i=0; i<W.rows(); i++){
      phi(i,0) = std::pow(t, i);
    }
    //std::cout << phi << "\n" << std::endl;
    return (a.transpose()*phi)(0,0);
  }
};
