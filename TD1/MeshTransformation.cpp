#include <cstdlib>
#include <igl/opengl/glfw/Viewer.h>
#include <ostream>

using namespace Eigen;

/**
 * A class for representing linear transformations on 3D points (using homogeneous coordinates)
 * */
class MeshTransformation{

    public:
        /**
         *Initialize the identity transformation
        */
        MeshTransformation(){
            MatrixXd m(4, 4);
            m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
            m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
            m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
            m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

            M = m;
        }


        /**
         *Initialize a scaling transformation
        */
        MeshTransformation(double s1, double s2, double s3){
            MatrixXd m(4,4);
            m << s1, 0, 0, 0,
                0, s2, 0, 0,
                0, 0, s3, 0,
                0, 0, 0, 1;
            M = m;
        }

        /**
         * Initialize a rotation transformation around a given axis (X, Y or Z) <br><br>
         * @param  direction  a value 0, 1 or 2 indicating the direction (X, Y or Z respectively)
        */
        MeshTransformation(double theta, int direction){
            MatrixXd m(4,4);
            switch (direction) {
                case 0: // around x axis
                    m << 1, 0, 0, 0,
                         0, cos(theta), -sin(theta), 0,
                         0, sin(theta), cos(theta), 0,
                         0, 0, 0, 1;
                    break;
                case 1: // around y axis
                    m << cos(theta), 0, sin(theta), 0,
                         0, 1, 0, 0,
                         -sin(theta), 0, cos(theta), 0,
                         0, 0, 0, 1;
                    break;
                case 2: // around z axis
                    m << cos(theta), -sin(theta), 0, 0,
                         sin(theta), cos(theta), 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 1;
                    break;
                default: // wrong input
                    exit(EXIT_FAILURE);
            }
            M = m;
        }

        /**
         * Initialize a translation
        */
        MeshTransformation(RowVector3d t) {
            MatrixXd m(4,4);
            m << 1, 0, 0, t.x(),
                0, 1, 0, t.y(),
                0, 0, 1, t.z(),
                0, 0, 0, 1;
            M = m;
        }

        /**
         * Matrix accessor
         * @return  the matrix transformation
        */
        MatrixXd get_matrix() {
            return M;
        }

        /**
         * Initialize a transformation given an input matrix 'm'
        */
        void set_matrix(MatrixXd m) {
            M = m;
        }

        /**
         * Apply the transformation to all vertices stored in a matrix 'V' <br>
         * @param V vector storing the input points
        */
        void transform(MatrixXd &V) {
            MatrixXd tmp = V;
            for(int i=0; i<V.rows(); i++){
                RowVector3d newRow = transform(V.row(i));
                tmp.row(i) = newRow;
            }
            V = tmp;
        }

        /**
         * Apply the transformation to a 3d (row) vector 'v' <br>
         * Remark: use homogeneous coordinates
         * @return  the vector after transformation
        */
        RowVector3d transform(RowVector3d v) {
            return (v.homogeneous()*M).hnormalized();
        }

        /**
        * Compose the current transformation with a transfomation 't': return a new transformation
        */
        MeshTransformation compose(MeshTransformation t) {
            MeshTransformation res(1.0, 1.0, 1.0);
            res.set_matrix(t.get_matrix()*M);
            return res;
        }

        /**
         * Print the matrix transformation
        */
        friend std::ostream& operator<<(std::ostream &os, MeshTransformation& t) {
            return os << "matrix:\n" << t.get_matrix() << std::endl;
        }

    private:
        MatrixXd M; // a 4x4 matrix representing a linear transformation

};
