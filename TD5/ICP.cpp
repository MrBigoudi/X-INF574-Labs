#include "ICP.h"
#include "octree.hpp"
#include <Eigen/src/Core/Matrix.h>

void brute_force_nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2){
    // return the nearest neighbour to V1 in V2 as nn_V2
    // Complete here
    // compute nearest neighbours
    nn_V2.resize(V1.rows(), V2.cols());
    for(int i=0; i<V1.rows(); i++){
        float minDist = -1.0f;
        for(int j=0; j<V2.rows(); j++){
            float curDist = (V1.row(i) - V2.row(j)).norm();
            if(minDist == -1 || curDist < minDist){
                nn_V2.row(i) = V2.row(j);
                minDist = curDist;
            }
        }
    }
}

void nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2){
    // return the nearest neighbour to V1 in V2 as nn_V2
    // Complete here
    // compute nearest neighbours
    Octree tree(V2);
    Eigen::MatrixXi indices = Eigen::MatrixXi::Zero(V2.rows(),1);
    tree.knn(V1, 1, indices);
    nn_V2.resize(V1.rows(), V2.cols());
    for(int i=0; i<V1.rows(); i++){
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
    nn_V2.resize(V1.rows(), V2.cols());

    // find the three closest point in V1 for each point in V1
    int k = 3;
    Octree tree(V1);
    Eigen::MatrixXi indices = Eigen::MatrixXi::Zero(V1.rows(), k);
    tree.knn(V1, k, indices);

    // for each point in V1
    for(int i=0; i<V1.rows(); i++){
        Eigen::Vector3d curPoint = V1.row(i);
        Eigen::Vector3d tmp1 = V1.row(indices(i, 1));
        Eigen::Vector3d tmp2 = V1.row(indices(i, 2));
        // build the plane
        Eigen::Vector3d v1 = tmp1 - curPoint;
        Eigen::Vector3d v2 = tmp2 - curPoint;
        Eigen::Vector3d normal = v1.cross(v2).normalized();
        // plane: ax + by + cz + d = 0
        float a = normal.x();
        float b = normal.y();
        float c = normal.z();
        float d = -normal.dot(curPoint);
        float deno = a * a + b * b + c * c;

        // for each point in V2, find the one closest to the plane
        float minDist = -1.0f;
        for(int j=0; j<V2.rows(); j++){
            // project the point in v2 on the plane
            Eigen::Vector3d p = V2.row(j);
            Eigen::Vector3d q = p - curPoint;
            float distAlongNorm = q.dot(normal);
            Eigen::Vector3d projectedPoint = p - distAlongNorm*normal;\
            // update the min distance and nn_V2
            float curDist = (curPoint - projectedPoint).norm();
            if(minDist < 0.0f || curDist < minDist){
                nn_V2.row(i) = p;
                minDist = curDist;
            }
        }
    }
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
