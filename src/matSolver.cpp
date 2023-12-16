#include "matSolver.h"
Mat_Solver::Mat_Solver()
    : PauliX((Eigen::MatrixXcd(2, 2) << std::complex<double>(0, 0), std::complex<double>(1, 0), std::complex<double>(1, 0), std::complex<double>(0, 0)).finished()),
      PauliY((Eigen::MatrixXcd(2, 2) << std::complex<double>(0, 0), std::complex<double>(0, -1), std::complex<double>(0, 1), std::complex<double>(0, 0)).finished()),
      PauliZ((Eigen::MatrixXcd(2, 2) << std::complex<double>(1, 0), std::complex<double>(0, 0), std::complex<double>(0, 0), std::complex<double>(-1, 0)).finished())
{
}
Eigen::VectorXcd Mat_Solver::basis(int n, int i)
{
    Eigen::VectorXcd basisVec(n);
    basisVec.setZero();
    basisVec[i] = std::complex<double>(1.0, 0.0); // 将第 i 个元素设置为复数 1.0 + 0.0i
    return basisVec;
}

std::vector<double> Mat_Solver::function(int n, double phi, double rho)
{
    Eigen::MatrixXcd result(n, n);
    result.setZero();
    auto start = std::chrono::high_resolution_clock::now();
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
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time for matrix computation: " << elapsed.count() << " s\n";
    auto start1 = std::chrono::high_resolution_clock::now();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(functionMat);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start1;
    std::cout << "Elapsed time for eigenvalue computation: " << elapsed1.count() << " s\n";
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

std::vector<double> Mat_Solver::function2(int n, double phi, double rho)
{
    Eigen::MatrixXcd result(2 * n, 2 * n);
    result.setZero();
    std::complex<double> img(0, 1);
// Parallel region for computing the first result matrix
#pragma omp parallel
    {
        Eigen::MatrixXcd local_result = Eigen::MatrixXcd::Zero(2 * n, 2 * n);

#pragma omp for nowait
        for (int i = 0; i < n - 1; i++)
        {
            // Eigen::MatrixXcd temp1 = Eigen::kroneckerProduct(basis(n, i + 1), basis(n, i).transpose());
            // Eigen::MatrixXcd temp2 = Eigen::kroneckerProduct(basis(n, i), basis(n, i + 1).transpose());
            local_result += 2 * (coeff1 * cos(M_PI * phi)) * Eigen::kroneckerProduct(PauliX, Eigen::kroneckerProduct(basis(n, i), basis(n, i).transpose())) +
                            0.5 * coeff2 * sin(M_PI * phi) * std::exp(rho * img) * Eigen::kroneckerProduct(PauliY, Eigen::kroneckerProduct(basis(n, i + 1), basis(n, i).transpose())) +
                            0.5 * coeff2 * sin(M_PI * phi) * std::exp(-rho * img) * Eigen::kroneckerProduct(PauliY, Eigen::kroneckerProduct(basis(n, i), basis(n, i + 1).transpose()));
        }

#pragma omp critical
        result += local_result;
    }

    result *= -img;
    result = result.exp();

    Eigen::MatrixXcd result2(2 * n, 2 * n);
    result2.setZero();

// Parallel region for computing the second result matrix
#pragma omp parallel
    {
        Eigen::MatrixXcd local_result2 = Eigen::MatrixXcd::Zero(2 * n, 2 * n);

#pragma omp for nowait
        for (int i = 0; i < n - 1; i++)
        {
            Eigen::MatrixXcd temp1 = Eigen::kroneckerProduct(basis(n, i), basis(n, i).transpose());
            Eigen::MatrixXcd temp2 = Eigen::kroneckerProduct(basis(n, i), basis(n, i + 1).transpose());
            Eigen::MatrixXcd temp3 = Eigen::kroneckerProduct(PauliZ, temp1);
            Eigen::MatrixXcd temp4 = Eigen::kroneckerProduct(PauliZ, temp2);
            local_result2 += (std::exp(-rho * img) * cos(M_PI * rho)) *
                                 Eigen::kroneckerProduct(PauliZ, Eigen::kroneckerProduct(basis(n, i), basis(n, i).transpose())) +
                             (0.5 * coeff3 * std::exp(rho * img) * cos(M_PI * rho)) * Eigen::kroneckerProduct(PauliZ, Eigen::kroneckerProduct(basis(n, i), basis(n, i + 1).transpose()));
        }

#pragma omp critical
        result2 += local_result2;
    }

    result2 *= -img;
    result2 = result2.exp();
    Eigen::MatrixXcd functionMat = result * result2;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(functionMat);

    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    std::vector<double> angles(eigenvalues.size());

    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        angles[i] = std::arg(eigenvalues[i]) / M_PI; // Store the phase angles in units of π
    }

    std::sort(angles.begin(), angles.end(), std::greater<double>());
    return angles;
}
