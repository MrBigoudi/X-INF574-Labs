#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F1; // incidence relations between faces and edges (f columns)

/**
 * This function is called every time a keyboard button is pressed
 * */
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1') {
    std::cout << "saving to OBJ format" << std::endl;
    igl::writeOBJ("../data/converted_mesh.obj", V1, F1);
  }

  if (key == '2')
  {
    MatrixXd transform(3,3);
    double zoomOutFactor = 1.2;
    transform << zoomOutFactor, 0, 0,
                  0, zoomOutFactor, 0,
                  0, 0, zoomOutFactor;

    for(int i=0; i<V1.rows(); i++){
      V1.row(i) = V1.row(i)*transform;
    }

    viewer.data(0).clear(); // Clear should be called before drawing the mesh
    viewer.data(0).set_mesh(V1, F1); // update the mesh (both coordinates and faces)
    //viewer.core().align_camera_center(V1, F1);
  }

  return false;
}

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V)
{
  // compute the corners of the bounding box
  Vector3d m = V.colwise().minCoeff();
  Vector3d M = V.colwise().maxCoeff();

  MatrixXd V_box(8, 3);  // Corners of the bounding box
  MatrixXi E_box(12, 2); // edges of the bounding box

  V_box << m(0), m(1), m(2),
      M(0), m(1), m(2),
      M(0), M(1), m(2),
      m(0), M(1), m(2),
      m(0), m(1), M(2),
      M(0), m(1), M(2),
      M(0), M(1), M(2),
      m(0), M(1), M(2);

  E_box << 0, 1,
      1, 2,
      2, 3,
      3, 0,
      4, 5,
      5, 6,
      6, 7,
      7, 4,
      0, 4,
      1, 5,
      2, 6,
      7, 3;

  viewer.append_mesh();
  viewer.data(1).add_points(V_box, Eigen::RowVector3d(1, 0, 0));

  for (unsigned i = 0; i < E_box.rows(); ++i) // Plot the edges of the bounding box
    viewer.data().add_edges(
        V_box.row(E_box(i, 0)),
        V_box.row(E_box(i, 1)),
        Eigen::RowVector3d(1, 0, 0));
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{
  igl::readOFF("../data/star.off", V1, F1); // Load an input mesh in OFF format

  //  print the number of mesh elements
  std::cout << "Vertices: " << V1.rows() << std::endl;
  std::cout << "Faces:    " << F1.rows() << std::endl;

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V1, F1); // load a face-based representation of the input 3d shape
  draw_bounding_box(viewer, V1); // draw the boundaing box (red edges and vertices)

  viewer.launch(); // run the editor
}
