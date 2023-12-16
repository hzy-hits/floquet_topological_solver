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
#include <chrono>

class Mat_Solver
{

public:
    Mat_Solver();
    void start();
    std::vector<double> function(int n, double phi, double rho);
    std::vector<double> function2(int n, double phi, double rho);

private:
    double coeff1 = 1;
    double coeff2 = 0.5;
    double coeff3 = 5;
    Eigen::VectorXcd basis(int n, int i);
    Eigen::MatrixXcd PauliX;
    Eigen::MatrixXcd PauliY;
    Eigen::MatrixXcd PauliZ;
};
