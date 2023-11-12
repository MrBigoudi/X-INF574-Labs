#ifndef __OCTREE_HPP__
#define __OCTREE_HPP__

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <igl/octree.h>
#include <igl/knn.h>

class Octree{
    private:
        std::vector<std::vector<int > > _PointIndices;
        Eigen::MatrixXi _Ch;
        Eigen::MatrixXd _Cn;
        Eigen::VectorXd _W;
        Eigen::MatrixXd _V;
    public:
        Octree(const Eigen::MatrixXd &V){
            _V = V;
            igl::octree(V, _PointIndices, _Ch, _Cn, _W);
        }

        void knn(const Eigen::MatrixXd &V, size_t k, Eigen::MatrixXi &I) const {
            igl::knn(V, _V, k, _PointIndices, _Ch, _Cn, _W, I);
        }
};

#endif