#include "matSolver.h"

Eigen::VectorXd Mat_Solver::basis(int n, int i)
{
    Eigen::VectorXd basisVec(n);
    basisVec.setZero();
    basisVec[i] = 1.0;
    return basisVec;
}

std::vector<double> Mat_Solver::function(int n, double phi, double rho)
{
    Eigen::MatrixXcd result(n, n);
    result.setZero();

// Parallel region for computing the first result matrix
#pragma omp parallel
    {
        Eigen::MatrixXcd local_result = Eigen::MatrixXcd::Zero(n, n);

#pragma omp for nowait
        for (int i = 0; i < n - 1; i++)
        {
            local_result += (coeff1 + pow(-1, i) * coeff2 * cos(M_PI * phi)) *
                                Eigen::kroneckerProduct(basis(n, i + 1), basis(n, i).transpose()) +
                            Eigen::kroneckerProduct(basis(n, i), basis(n, i + 1).transpose());
        }

#pragma omp critical
        result += local_result;
    }
    std::complex<double> i(0, 1);
    result *= -i;
    result = result.exp();

    Eigen::MatrixXcd result2(n, n);
    result2.setZero();

// Parallel region for computing the second result matrix
#pragma omp parallel
    {
        Eigen::MatrixXcd local_result2 = Eigen::MatrixXcd::Zero(n, n);

#pragma omp for nowait
        for (int i = 0; i < n; i++)
        {
            local_result2 += (pow(-1, i) * coeff3 * cos(M_PI * rho)) *
                             Eigen::kroneckerProduct(basis(n, i), basis(n, i).transpose());
        }

#pragma omp critical
        result2 += local_result2;
    }

    result2 *= -i;
    result2 = result2.exp();

    Eigen::MatrixXcd functionMat = result * result2;

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(functionMat);
    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    std::vector<double> angles(eigenvalues.size());

    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        angles[i] = std::arg(eigenvalues[i]) / M_PI; // Store the phase angles in units of π
    }

    // Sort the phase angles
    std::sort(angles.begin(), angles.end(), std::greater<double>());

    // Output the sorted phase angles
    // std::cout << "Sorted angles of eigenvalues (in π units):" << std::endl;
    // for (const auto &angle : angles)
    // {
        // std::cout << angle << std::endl;
    // }
    // std::cout << angles.size() << std::endl;

    return angles;
}
