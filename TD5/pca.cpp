#include "pca.h"
#include "octree.hpp"


void k_nearest_neighbour(const MatrixXd &V1,Eigen::MatrixXi &I, int k){
    // return the k nearest neighbour index
    // complete here
    Octree tree(V1);
    tree.knn(V1, k, I);
}



void compute_normals(const MatrixXd &V1,const Eigen::MatrixXi &I, int k, MatrixXd &normals){
    // compute the normals using PCA
    normals = Eigen::MatrixXd::Zero(V1.rows(), 3);

    // for each point, compute its covariance matrix
    for(int p=0; p<V1.rows(); p++){
        // get the current point
        Eigen::MatrixXd curPoint = V1.row(p);

        // get its neighbours
        Eigen::MatrixXd neighbours = Eigen::MatrixXd::Zero(k,3);
        for(int i=0; i<k; i++){
            neighbours.row(i) = V1.row(I(p, i));
        }

        // get the mean of the neighbours
        Eigen::MatrixXd pMean = neighbours.colwise().mean();
        
        //build the covarianve matrix
        Eigen::MatrixXd crossCovariance = Eigen::MatrixXd::Zero(3,3);
        for(int i=0; i<neighbours.rows(); i++){
            crossCovariance += (neighbours.row(i)-pMean).transpose()*(neighbours.row(i)-pMean);
        }

        // find the first eigenvector
        Eigen::EigenSolver<Eigen::MatrixXd> solver(crossCovariance);
        Eigen::VectorXcd eigenvalues = solver.eigenvalues();
        int smallestIndex = 0;
        for (int i = 1; i < eigenvalues.size(); i++) {
            if (eigenvalues[i].real() < eigenvalues[smallestIndex].real()) {
                smallestIndex = i;
            }
        }
        Eigen::Vector3d normal = solver.eigenvectors().col(smallestIndex).real().cast<double>();
        normals.row(p) = normal.normalized() /  5.0f; // divide to 
    }

}
