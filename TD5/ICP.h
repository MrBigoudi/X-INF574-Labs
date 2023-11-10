
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
using namespace Eigen;
#include <igl/octree.h>
#include <igl/knn.h>

void nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2);
void nearest_neighbour_point_to_plane(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2);
void transform(MatrixXd &V1,const MatrixXd &V2);
float getSumPairwiseNN(const MatrixXd& V1, const MatrixXd& V2);