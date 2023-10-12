#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include "MeshTransformation.cpp"
#include <iostream>
#include <ostream>

using namespace Eigen;

MatrixXd V1;
MatrixXi F1;

MeshTransformation *identity = new MeshTransformation();         // identity transformation
MeshTransformation scalingUp(1.2, 1.2, 1.2);                     // scaling (increase size)
MeshTransformation scalingDown(0.9, 0.9, 0.9);                   // scaling (decrease size)
MeshTransformation translationPlusZ(RowVector3d(0., 0., 0.2));   // translation in z-direction (positive)
MeshTransformation translationMinusZ(RowVector3d(0., 0., -0.2)); // translation in z-direction (negative)
MeshTransformation rotationX(0.2, 0);                            // rotation around X axis

MeshTransformation *T = identity; // by default set the identity transformation

/**
 * Compute the barycenter of 'n' points
 * @param V  a matrix (storing double values) having 'n' rows and '3' columns
 * @return  a row vector (3 columns) corresponding to the barycenter of the input points
 * */
RowVector3d compute_barycenter(MatrixXd &V){
    RowVector3d bar(0., 0., 0.);
    int nbPoints = V.rows();
    //std::cout << "nbPoints: " << nbPoints << std::endl;
    if(nbPoints == 0) return bar; // if no points given

    for(int i=0; i<nbPoints; i++){
        RowVector3d rowTmp = V.row(i);
        //std::cout << "rowTmp: " << rowTmp << std::endl;
        bar += rowTmp;
    }
    //std::cout << bar << std::endl;

    return bar / nbPoints;
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier){
    std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
    if (key == 'B')
    {
        T = &scalingUp;
    }
    else if (key == 'S')
    {
        T = &scalingDown;
    }
    else if (key == 'U')
    {
        T = &translationPlusZ;
    }
    else if (key == 'D')
    {
        T = &translationMinusZ;
    }
    else if (key == 'X')
    {
        T = &rotationX;
    }
    else
        return false;

    // apply the transformation T to all points in V1
    // Warning: the rotation center should coincide with the barycenter of V1
    
    // move the barycenter to the origin
    // do the transformations
    // move the barycenter back in its place
    MeshTransformation centerBar(-compute_barycenter(V1));
    MeshTransformation restoreBar(compute_barycenter(V1));
    centerBar.compose(T->compose(restoreBar)).transform(V1);

    viewer.data(0).clear();
    viewer.data(0).set_mesh(V1, F1);

    return false;
}

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V){
    // compute the corners of the bounding box
    Vector3d m = V.colwise().minCoeff();
    Vector3d M = V.colwise().maxCoeff();

    MatrixXd V_box(8, 3);  // Corners of the bounding box
    MatrixXi E_box(12, 2); // edges of the boundaing box

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
int main(int argc, char *argv[]){
    igl::readOFF("../data/star.off", V1, F1);

    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = &key_down;
    viewer.data().set_mesh(V1, F1);

    draw_bounding_box(viewer, V1);
    viewer.launch();
}