#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;
class HermiteInterpolation{

  MatrixXd solutionx;// solution for x
  MatrixXd solutiony;// solution for y
  MatrixXd steps;// points to interpolate
  MatrixXd solvey;// left size of linear system
  MatrixXd solvex;// left size of linear system
  MatrixXd slopes;
  MatrixXd W;

void solve_x(){
  MatrixXd inv = W.inverse();
  // for every functions fi
  for(int i=0; i<steps.rows()-1; i++){
    // create the current solutions
    solutionx(0,0) = steps.row(i).x();
    solutionx(1,0) = steps.row(i+1).x();
    solutionx(2,0) = slopes.row(i).x();
    solutionx(3,0) = slopes.row(i+1).x();
    // get the temporary solutions
    MatrixXd tmp = inv*solutionx;
    // fill the solvex matrix
    int idx = i*4;
    solvex(idx,0)   = tmp(0,0);
    solvex(idx+1,0) = tmp(1,0);
    solvex(idx+2,0) = tmp(2,0);
    solvex(idx+3,0) = tmp(3,0);
  }
}

void solve_y(){
  MatrixXd inv = W.inverse();
  // for every functions fi
  for(int i=0; i<steps.rows()-1; i++){
    // create the current solutions
    solutiony(0,0) = steps.row(i).y();
    solutiony(1,0) = steps.row(i+1).y();
    solutiony(2,0) = slopes.row(i).y();
    solutiony(3,0) = slopes.row(i+1).y();
    // get the temporary solutions
    MatrixXd tmp = inv*solutiony;
    // fill the solvex matrix
    int idx = i*4;
    solvey(idx,0)   = tmp(0,0);
    solvey(idx+1,0) = tmp(1,0);
    solvey(idx+2,0) = tmp(2,0);
    solvey(idx+3,0) = tmp(3,0);
  }
}

public:
  HermiteInterpolation(const MatrixXd &V1){
    W = MatrixXd::Zero(4,4);
    // complete here
    W << 
        1, 0, 0, 0,
        1, 1, 1, 1,
        0, 1, 0, 0,
        0, 1, 2, 3;

    steps = V1;

    slopes = MatrixXd::Zero(5,2);
    slopes(0,0) = 1; slopes(0,1) = 1;
    slopes(1,0) = 1; slopes(1,1) = 1;
    slopes(2,0) = 1.0; slopes(2,1) = 1;
    slopes(3,0) = -2.5; slopes(3,1) = 1;
    slopes(4,0) = 1.0; slopes(4,1) = 1;

    // complete here
    solutionx = MatrixXd::Zero(4,1);
    solutiony = MatrixXd::Zero(4,1);
    solvex = MatrixXd::Zero(4*(steps.rows()-1),1);
    solvey = MatrixXd::Zero(4*(steps.rows()-1),1);

    solve_x();
    solve_y();

    // std::cout
    //   << "\nsteps:\n" << steps
    //   << "\nslopes:\n" << slopes
    //   << "\nW:\n" << W
    //   << "\nsolutionx:\n" << solutionx
    //   << "\nsolutiony:\n" << solutiony
    //   << "\nsolvex:\n" << solvex
    //   << "\nsolvey:\n" << solvey
    //   << std::endl;
  }

  // complete linspace with coordinates
  void eval_function(MatrixXd &linspace){
    // complete here
    int resolution = linspace.rows();
    int step = resolution / (steps.rows()-1);

    int cpt = 0;
    float beg = 0.0f;
    float end = 1.0f;
    float dt = ((end-beg)/step);

    // for every piecewise curve
    for(int curve = 0; curve < steps.rows()-1; curve++){
      int idx = curve*4;
      float a_xi = solvex(idx,0); 
      float b_xi = solvex(idx+1,0); 
      float c_xi = solvex(idx+2,0); 
      float d_xi = solvex(idx+3,0); 

      float a_yi = solvey(idx,0); 
      float b_yi = solvey(idx+1,0); 
      float c_yi = solvey(idx+2,0); 
      float d_yi = solvey(idx+3,0); 

      // create step points between these two values
      for(float t=beg; t<end-dt; t+=dt){
        float x = a_xi + b_xi*t + c_xi*t*t + d_xi*t*t*t;
        float y = a_yi + b_yi*t + c_yi*t*t + d_yi*t*t*t;
        // std::cout << "t: " << t 
        //           << ", cpt: " << cpt 
        //           << ", x: " << x 
        //           << ", y: " << y 
        //           << std::endl;

        linspace(cpt, 0) = x;
        linspace(cpt, 1) = y;
        cpt++;
      }
    }
  }

  /*
  Evaluate tangent at step i
  */
  void eval_tangent(float i, MatrixXd &dX, float t){
    // complete here
    float dt = 0.25f;
    // calculate fxi
    float a_xi = solvex(4*i,0);
    float b_xi = solvex(4*i+1,0);
    float c_xi = solvex(4*i+2,0);
    float d_xi = solvex(4*i+3,0);
    // calculate fyi
    float a_yi = solvey(4*i,0);
    float b_yi = solvey(4*i+1,0);
    float c_yi = solvey(4*i+2,0);
    float d_yi = solvey(4*i+3,0);

    // y = f(a) + f'(a)(x-a)
    float fx = a_xi + b_xi*t + c_xi*t*t + d_xi*t*t*t;
    float fy = a_yi + b_yi*t + c_yi*t*t + d_yi*t*t*t;
    float fx_prime = b_xi + 2*c_xi*t + 3*d_xi*t*t;
    float fy_prime = b_yi + 2*c_yi*t + 3*d_yi*t*t;
    dX(i, 0) = fx + fx_prime*dt;
    dX(i, 1) = fy + fy_prime*dt;
  }
};
