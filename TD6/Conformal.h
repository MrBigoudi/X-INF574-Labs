#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

class ConformalParametrization {

    private:
        VectorXi boundary_flat;
    public:

        MatrixXd V_uv,V_uv_harmonic, V,bc;
		MatrixXi F;
		SparseMatrix<double> Dirichlet, Area, ConformalEnergy;
		
		VectorXi boundary,b;
        ConformalParametrization(const MatrixXd &, const MatrixXi &);

        double cotangent(double );
        //compute the matrix for dirichlet energy of the system
        void compute_dirichlet(SparseMatrix<double> &Dirichlet);

        void compute_harmonic_parametrization(MatrixXd &V_uv_harmonic);

        void compute_area(SparseMatrix<double> &Area);

        void compute_conformal_energy(SparseMatrix<double> &ConformalEnergy);

        void minimize_energy(MatrixXd &V_uv);

        void build_parametrizations();

        void compute_conformal_factor(MatrixXd&);

        void compute_conformal_coloring(const MatrixXd &, MatrixXd );
};

