#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class LinearInterpolation{
  MatrixXd V;
  public:
    LinearInterpolation(const MatrixXd &V0){
      V = V0;
    }

    // evaluate function at time t
    float eval_function(float t){
      // get the point left to t
      RowVector3d leftPoint;
      RowVector3d rightPoint;
      //std::cout << leftPoint << ", " << rightPoint << "\n" << std::endl;
      for(int i=0; i<V.rows(); i++){
        if(t<V.row(i).x()){
          rightPoint = V.row(i);
          leftPoint = V.row(i-1);
          break;
        }
      }

      float a = (leftPoint.y()-rightPoint.y()) / (leftPoint.x()-rightPoint.x());
      float b = leftPoint.y() - a*leftPoint.x();

      return a*t + b;
    }
  };
