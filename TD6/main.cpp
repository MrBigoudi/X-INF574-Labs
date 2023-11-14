#include <igl/boundary_loop.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/lscm.h>
#include "Conformal.h"



Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd V_uv_libigl;
Eigen::MatrixXd V_uv_your;
Eigen::MatrixXd V_uv_harmonic;



bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{

  if (key == '1')
  {
    // Plot the 3D mesh
    viewer.data().set_mesh(V,F);
    viewer.core().align_camera_center(V,F);
  }
  else if (key == '2')
  {
    // Plot the mesh in 2D using the UV coordinates as vertex coordinates
    viewer.data().set_mesh(V_uv_your, F);
    viewer.data().set_uv(V_uv_your);
    viewer.core().align_camera_center(V_uv_your, F);
    
  }
  else if (key == '3') {
      viewer.data().set_mesh(V_uv_harmonic, F);
      viewer.data().set_uv(V_uv_harmonic);
      viewer.core().align_camera_center(V_uv_harmonic, F);
  }
  else if (key == '4') {
      viewer.data().set_mesh(V_uv_libigl, F);
      viewer.data().set_uv(V_uv_libigl);
      viewer.core().align_camera_center(V_uv_libigl, F);
  }
  

  viewer.data().compute_normals();

  return false;
}



int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OFF format
  
  //use this line for Visual Studio Code and the CMake PlugIn
  // igl::readOFF("../../data/camelhead.off", V, F);

  //use this line for Visual Studio
  //igl::readOFF("../../../data/beetle.off", V, F);

  //use this line for Linux or Mac
  igl::readOFF("./data/camelhead.off", V, F);
  // igl::readOFF("../data/beetle.off", V, F);

  VectorXi bnd,b(2,1);
  igl::boundary_loop(F,bnd);
  b(0) = bnd(0);
  b(1) = bnd(bnd.size()/2);
  MatrixXd bc(2,2);
  bc<<0,0,1,0;
  //use the in-build Libigl parametrization to compare  
  igl::lscm(V,F,b,bc,V_uv_libigl);

  // Scale the uv
  V_uv_libigl *= 5;

  //Now set up your own system
  ConformalParametrization C_sys(V,F);

  C_sys.build_parametrizations();
  V_uv_harmonic = C_sys.V_uv_harmonic;
  V_uv_your = C_sys.V_uv;

  V_uv_your*=5;

  //--------Uncomment this part for the coloring of the conformal factor
  // MatrixXd conformal_factor;
  // MatrixXd Color;
  // C_sys.compute_conformal_factor(conformal_factor);
  // C_sys.compute_conformal_coloring(conformal_factor,Color);
  //Now: Set the color on the mesh!
  

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_uv(V_uv_libigl);
  viewer.callback_key_down = &key_down;

  // Disable wireframe
  viewer.data().show_lines = false;

  // Draw checkerboard texture
  viewer.data().show_texture = true;

  // Launch the viewer
  viewer.launch();
}