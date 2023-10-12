#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class CubicInterpolation{
  MatrixXd M;// coefficients
  MatrixXd a;// solution
  MatrixXd V;// points to interpolate
  MatrixXd y;// left size of linear system

  /*
  Initialize system constraints
  */
  void init_system(const MatrixXd &V, MatrixXd &M, MatrixXd &y){
    // complete here
    int rows = 4*(V.rows()-1);
    int i = 0;
    int cpt = 0;
    for(int j=0; j<=rows; j+=4){
      float yi = V.row(cpt).y();
      float xi = V.row(cpt).x();
      //std::cout << "cpt: " << cpt << ", j: " << j << ", yi: " << yi << std::endl;

      // create the ys
      if(j==0){
        y(j,0) = 0;
        y(j+1,0) = yi;
      }
      else if(j==rows) {
        y(j-2,0) = 0;
        y(j-2+1,0) = yi;
      }
      else{
        y(j-2,0) = yi;
        y(j-2+1,0) = 0;
        y(j-2+2,0) = 0;
        y(j-2+3,0) = yi;
      }

      // ceate the Ms
      int ai_1 = i-4, bi_1 = i-4+1, ci_1 = i-4+2, di_1 = i-4+3;
      int ai = i, bi = i+1, ci = i+2, di = i+3;

      // if i == 0, only f(x0) = y0 and f'(x0) = 0
      // if i == 4(n-1)-2, only f(xn-1) = yn-1 and f'(xn-1) = 0
      if(j==0){
        M(j, bi) = 1;
        M(j, ci) = 2*xi;
        M(j, di) = 3*xi*xi;
        M(j+1, ai) = 1;
        M(j+1, bi) = xi;
        M(j+1, ci) = xi*xi;
        M(j+1, di) = xi*xi*xi;
        i+=4;
        cpt++;
        continue;
      } else if(j==rows){
        //f'n-1(xn-1) = 0
        M(j-2, bi_1) = 1;
        M(j-2, ci_1) = 2*xi;
        M(j-2, di_1) = 3*xi*xi;
        //fn-1(xn-1) = yn-1
        M(j-2+1, ai_1) = 1;
        M(j-2+1, bi_1) = xi;
        M(j-2+1, ci_1) = xi*xi;
        M(j-2+1, di_1) = xi*xi*xi;
        break;
      }

      // fi-1(xi) = yi
      // ai-1 + bi-1xi + ci-1xi^2 + di-1xi^3
      M(j-2, ai_1) = 1; 
      M(j-2, bi_1) = xi; 
      M(j-2, ci_1) = xi*xi;
      M(j-2, di_1) = xi*xi*xi;

      // f'i-1(xi) - f'i(xi) = 0
      // bi-1 + 2ci-1xi + 3di-1xi^2 - bi - 2cixi - 3dixi^2
      M(j-2+1, bi_1) = 1;
      M(j-2+1, ci_1) = 2*xi;
      M(j-2+1, di_1) = 3*xi*xi;
      M(j-2+1, bi) = -1;
      M(j-2+1, ci) = -2*xi;
      M(j-2+1, di) = -3*xi*xi;

      // f"i-1(xi) - f"1(xi) = 0
      // 2ci-1 + 6di-1xi - 2ci - 6dixi
      M(j-2+2, ci_1) = 2;
      M(j-2+2, di_1) = 6*xi;
      M(j-2+2, ci) = -2;
      M(j-2+2, di) = -6*xi;

      // fi(xi) = yi
      // ai + bixi+1 + cixi+1^2 + dixi+1^3
      M(j-2+3, ai) = 1; 
      M(j-2+3, bi) = xi; 
      M(j-2+3, ci) = xi*xi;
      M(j-2+3, di) = xi*xi*xi;

      i += 4;
      cpt++;
    }
  }


public:
  CubicInterpolation(const MatrixXd &V1){
    M = MatrixXd::Zero(4*(V1.rows() - 1),4*(V1.rows() - 1));
    y = MatrixXd::Zero(4*(V1.rows() - 1),1);
    V = V1;
    //std::cout << "before init V:\n" << V << "\n" << std::endl;
    init_system(V1, M, y );
    // complete here
    //std::cout << "M:\n" << M << "\n" << std::endl;
    //std::cout << "y:\n" << y << "\n" << std::endl;
    a = M.inverse()*y;
    //std::cout << "a:\n" << a << "\n" << std::endl;
  }


  /*
  Evaluate tangent at step i
  */
  void eval_tangent(float i, MatrixXd &dX, float x){
    // complete here
    float dx = 0.25f;
    // calculate fi(xi)
    float ai = a(4*i,0);
    float bi = a(4*i+1,0);
    float ci = a(4*i+2,0);
    float di = a(4*i+3,0);

    // y = f(a) + f'(a)(x-a)
    float f = ai + bi*x + ci*x*x + di*x*x*x;
    float f_prime = bi + 2*ci*x + 3*di*x*x;
    dX(i, 0) = x+dx;
    dX(i, 1) = f + f_prime*dx;
  }

  /*
  Evaluate function at time t
  */
  float eval_function(float t){
    // complete here
    // get the point left to t
    int i = 0;
    bool found = false;
    int rows = V.rows();
    //std::cout << "final i: " << i << ", final i+1" << i_plus_1 << "\n" << std::endl;
    for(int j=0; j<rows; j++){
      if(t<V.row(j).x()){
        i = j-1;
        found = true;
        break;
      }
    }
    // calculate fi(xi)
    if(!found) return y(rows-2+1,0);
    float ai = a(4*i,0);
    float bi = a(4*i+1,0);
    float ci = a(4*i+2,0);
    float di = a(4*i+3,0);
    return ai + bi*t + ci*t*t + di*t*t*t;
  }
};
