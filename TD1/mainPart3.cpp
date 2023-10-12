#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <iostream>
#include <ostream>

using namespace Eigen;

MatrixXd V1; // vertex coordinates of the input mesh
MatrixXi F1; // incidence relations between faces and edges
RowVector3d rotation_axis(1., 0., 0.); // rotation axis

  	/**
	 * Print the components of a quaternion
	 */
	void print_quaternion(Quaterniond q) {
    std::cout << "(" << q.w() << ", " << q.x() << ", "<< q.y() << ", "<< q.z() << "\n";
	}

  	/**
	 * Apply to the input points (stored in a matrix V) a rotation of angle 'theta' around a given direction defined by vector 'u'
   * 
   * @param V  a matrix storing the input points ('n' rows, '3' columns)
   * @param u  a vector corresponding to the rotation axis
   * @param theta  rotation angle
	 */
  void transform(MatrixXd &V, RowVector3d u, double theta) {
    // for some reasons, u.normalize() is not working
    double sumU = u.x()+u.y()+u.z();
    RowVector3d uUnit(u.x()/sumU, u.y()/sumU, u.z()/sumU);
    Quaterniond q(cos(theta/2.0), sin(theta/2)*uUnit.x(), sin(theta/2)*uUnit.y(), sin(theta/2)*uUnit.z());
    Quaterniond qStar = q.inverse();
    int nbPoints = V.rows();
    for(int i=0; i<nbPoints; i++){
        Quaterniond v(0, V.row(i).x(), V.row(i).y(), V.row(i).z());
        v = q*v*qStar;
        V.row(i)[0] = v.x();
        V.row(i)[1] = v.y();
        V.row(i)[2] = v.z();
    }
  }

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;

    if (key == 'Q')
    {
      std::cout << "rotating with quaternions" << std::endl;
      MatrixXd &V = V1;
      MatrixXi &F = F1;
      RowVector3d u(1., 0., 0.); // rotation axis

      double angle=0.1;
      transform(V, u, angle); // apply the rotation around vector 'u' to all points in V

      viewer.data(0).clear();
      viewer.data(0).set_mesh(V, F);
      //viewer.core(mesh_selection).align_camera_center(V, F);
    }

    if (key == 'R') // generate a random direction
    {
      rotation_axis(0)=random(); // x component of the rotation axis
      rotation_axis(1)=random(); // y component
      rotation_axis(2)=random(); // z component
    }

  return false;
}

// This function is called every time a keyboard button is pressed
bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
      //std::cout << "rotating with quaternions" << std::endl;
      MatrixXd &V = V1;
      MatrixXi &F = F1;

      transform(V, rotation_axis, 0.05); // rotate the object

      viewer.data(0).clear();
      viewer.data(0).set_mesh(V, F);

  return false;
}

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V) {
  // compute the corners of the bounding box
  Vector3d m = V.colwise().minCoeff();
  Vector3d M = V.colwise().maxCoeff();

  MatrixXd V_box(8,3); // Corners of the bounding box
  MatrixXi E_box(12,2); // edges of the bounding box

  V_box <<
  m(0), m(1), m(2),
  M(0), m(1), m(2),
  M(0), M(1), m(2),
  m(0), M(1), m(2),
  m(0), m(1), M(2),
  M(0), m(1), M(2),
  M(0), M(1), M(2),
  m(0), M(1), M(2);
  
  E_box <<
  0, 1,
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
  7 ,3;

  viewer.append_mesh();
  viewer.data(1).add_points(V_box,Eigen::RowVector3d(1,0,0));
    
  for (unsigned i=0;i<E_box.rows(); ++i) // Plot the edges of the bounding box
    viewer.data().add_edges
    (
      V_box.row(E_box(i,0)),
      V_box.row(E_box(i,1)),
      Eigen::RowVector3d(1,0,0)
    );

}

// ------------ main program ----------------
int main(int argc, char *argv[]) {
  igl::readOFF("../data/star.off", V1, F1); // Load an input mesh in OFF format

  // input mesh
  std::cout << "Vertices: " << V1.rows() << std::endl;
  std::cout << "Faces:    " << F1.rows() << std::endl;

  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.callback_pre_draw = &pre_draw;
  viewer.data().set_mesh(V1, F1);
  viewer.core().is_animating = true; // animation is active

  draw_bounding_box(viewer, V1);
  viewer.launch();
}