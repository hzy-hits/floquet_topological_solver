#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <complex>


class Mat_Solver
{
public:
    void start();
    std::vector<double> function(int n, double phi, double rho);

private:
    double coeff1 = 1;
    double coeff2 = 0.5;
    double coeff3 = 5;
    Eigen::VectorXd basis(int n, int i);
};
