#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include "lagrange.cpp"
#include "linear.cpp"
#include "cubic.cpp"
#include "hermite.cpp"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;

/* Draw a curve given a list of points. Successive points in the list are connected by a segment*/
void draw_curve(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V){
  viewer.append_mesh();
  for (unsigned i = 0; i < V.rows()-1; ++i)
    viewer.data().add_edges(
        V.row(i),
        V.row(i+1),
        Eigen::RowVector3d(1, 0, 0));
}

/* Draw tangents on curve */
void draw_tangents(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V, const MatrixXd &dX){
  // complete here
  viewer.append_mesh();
  for (unsigned i = 1; i < V.rows()-1; i++){
    viewer.data().add_edges(
        Eigen::RowVector3d(V.row(i).x(), V.row(i).y(), 0),
        Eigen::RowVector3d(dX(i,0), dX(i, 1), 0),
        Eigen::RowVector3d(0, 0, 1));
  }
}

/*draw points from the list of points V*/
void draw_points(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V){
  viewer.append_mesh();
  viewer.data(0).add_points(V, Eigen::RowVector3d(1, 0, 0));
}

/*interpolate the leftmost and rightmost X position*/
void build_linspace(MatrixXd &linspace,const MatrixXd &V){
  for (size_t i = 0; i < linspace.rows(); i++) {
    linspace(i, 0) = V.col(0).minCoeff() + ((V.col(0).maxCoeff() - V.col(0).minCoeff())/(linspace.rows() - 1))*i;
  }
}

/*
int main(int argc, char *argv[])
{
  
  //  use this line for linux and MAC   
  igl::readPLY("../data/curve0.ply", V1, F1);
  
  //  use this line for visual sudio
  //  igl::readPLY("../../../data/curve0.ply", V1, F1);
  //  print the number of mesh elements
  std::cout << "Points: " << V1.rows() << std::endl;

  // choose interpolation method
  //LinearInterpolation interp(V1);
  //LagrangeInterpolation interp(V1);
  CubicInterpolation interp(V1);

  int resolution = 500; // resolution of displayed curve
  MatrixXd linspace = MatrixXd::Zero(resolution, 3);
  build_linspace(linspace, V1); // initialize the X axis of the interpolation

  for (size_t i = 0; i < resolution; i++) {
    linspace(i, 1) = interp.eval_function(linspace(i, 0));
  }

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  draw_points(viewer, V1); // draw the bounding box (red edges and vertices)
  draw_curve(viewer, linspace);

  // uncomment if we want to draw tangents
  MatrixXd dX(V1.rows(), 2);
  for (size_t i = 1; i < V1.rows()-1; i++) {
    interp.eval_tangent(i, dX, V1(i, 0));
  }
  draw_tangents(viewer, V1, dX);
  viewer.launch(); // run the editor
}
*/



// use this main for exercice 3
int main(int argc, char *argv[]){
  igl::readPLY("../data/curve0.ply", V1, F1);
  std::cout << "Points: " << V1.rows() << std::endl;
  HermiteInterpolation interp(V1);
  int resolution = 400;
  MatrixXd linspace = MatrixXd::Zero(resolution, 3);
  interp.eval_function(linspace); // eval function updates the entirety of linspace
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  draw_points(viewer, V1); // draw the bounding box (red edges and vertices)
  draw_curve(viewer, linspace);

  // uncomment if we want to draw tangents
  MatrixXd dX(V1.rows(), 2);
  for (size_t i = 1; i < V1.rows()-1; i++) {
    interp.eval_tangent(i, dX, 0.f);
  }
  draw_tangents(viewer, V1, dX);
  viewer.launch(); // run the editor
  return 0;
}