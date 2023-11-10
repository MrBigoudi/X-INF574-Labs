#include "ICP.h"
#include <Eigen/src/Core/Matrix.h>

class Octree{
    private:
        std::vector<std::vector<int > > _PointIndices;
        Eigen::MatrixXi _Ch;
        Eigen::MatrixXd _Cn;
        Eigen::VectorXd _W;
        Eigen::MatrixXd _V;
    public:
        Octree(const MatrixXd &V){
            _V = V;
            igl::octree(V, _PointIndices, _Ch, _Cn, _W);
        }

        void knn(const MatrixXd &V, size_t k, Eigen::VectorXi &I) const {
            igl::knn(V, _V, k, _PointIndices, _Ch, _Cn, _W, I);
        }
};

void nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2){
    // return the nearest neighbour to V1 in V2 as nn_V2
    // Complete here
    // compute nearest neighbours
    Octree tree(V1);
    Eigen::VectorXi indices = Eigen::VectorXi::Zero(V2.rows(),1);
    tree.knn(V2, 1, indices);
    nn_V2.resize(V2.rows(), V2.cols());
    for(int i=0; i<V2.rows(); i++){
        nn_V2.row(i) = V2.row(indices(i,0));
    }
}

float getSumPairwiseNN(const MatrixXd& V1, const MatrixXd& V2){
    float sum = 0.0f;
    for(int i=0; i<V1.rows(); i++){
        sum += (V1.row(i) - V2.row(i)).norm();
    }
    return sum;
}

void nearest_neighbour_point_to_plane(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2){
    // return the nearest neighbour to V1 in V2 as nn_V2 using the point to plane algorithm
    // Complete here
}




void transform(MatrixXd &V1,const MatrixXd &V2){
    //align V1 to V2 when V1 and V2 points are in correspondance
    // complete here

    // do translation
    Eigen::MatrixXd v1Mean = V1.colwise().mean();
    Eigen::MatrixXd v2Mean = V2.colwise().mean();
    Eigen::MatrixXd tmp = v2Mean-v1Mean;
    Eigen::MatrixXd translation(4,4);
    translation << 1., 0., 0., tmp(0,0),
                   0., 1., 0., tmp(0,1),
                   0., 0., 1., tmp(0,2),
                   0., 0., 0., 1.0;

    // build cross covariance matrix
    Eigen::MatrixXd crossCovariance = Eigen::MatrixXd::Zero(3,3);
    for(int i=0; i<V1.rows(); i++){
        crossCovariance += (V1.row(i)-v1Mean).transpose()*(V2.row(i)-v2Mean);
    }
    
    // svd
    JacobiSVD<MatrixXd> svd(crossCovariance, ComputeThinU | ComputeThinV);
    MatrixXd u = svd.matrixU();
    MatrixXd v = svd.matrixV();

    // do rotation
    Eigen::MatrixXd rotationTmp = v*u.transpose();
    if(rotationTmp.determinant() <= 0){
        v.col(2) *= -1;
        rotationTmp = v*u.transpose();
    }
    Eigen::MatrixXd rotation = Eigen::MatrixXd::Zero(4, 4);
    rotation.topLeftCorner(3, 3) = rotationTmp;
    rotation.row(3) << 0., 0., 0., 1.;

    // update V1
    for(int i=0; i<V1.rows(); i++){
        tmp = Eigen::MatrixXd::Zero(4,1);
        tmp.topLeftCorner(3, 1) = V1.row(i).transpose();
        tmp.row(3) << 1.;
        tmp = rotation*translation*tmp;
        V1.row(i) = tmp.topLeftCorner(3,1).transpose();
    }
}
