#include <iostream>

#include "matSolver.h"
#include <fstream>
#include <string>
int main(int, char **)

{
    std::cout << "Hello, from Mat_Solver!\n";

    Mat_Solver solver;
    std::string filename = "../data/data.txt";
    std::ofstream file;
    file.open(filename);

    std::vector<double> result;
    for (float phi1 = -1.1; phi1 <= 1.1; phi1 += 0.1)
    {
        for (float rho1 = -1.1; rho1 <= 1.1; rho1 += 0.1)

        {
            auto start = std::chrono::high_resolution_clock::now();
            file << phi1 << " " << rho1 << " " << std::endl;
            result = solver.function2(100, phi1, rho1);
            for (int i = 0; i < result.size(); i++)
            {
                if (file.is_open())
                {
                    file << result[i] << " ";
                }
                else
                {
                    std::cout << "Unable to open file";
                    return -1;
                }
            }
            file << std::endl;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        }
    }
    file.close();
}
