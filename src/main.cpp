#include <iostream>
#include "matSolver.h"

int main(int, char **)

{
    std::cout << "Hello, from Mat_Solver!\n";
    Mat_Solver solver;
    std::vector<double> result;
    result = solver.function(10, 1, 1);
}
