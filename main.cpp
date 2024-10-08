#include "utils/Vector.hpp"
// #include "solvers/JacobiSequential.hpp"
#include "solvers/jacobi_MPI.hpp"
#include "solvers/jacobi_sequential.hpp"
// #include "solvers/jacobi_seq.hpp"



// double a = 20.;
// double b = 2.;
// double U_0 = 0.1;
// double alpha = 0.5;

// double k_y = 2. * M_PI / b;
// double k_x = M_PI / (2. * a);
// double V(double y)
// {
//     return 1 - cos(k_y * y);
// }

// double V_tilde(double x)
// {
//     return 1 - sin(k_x * x);
// }

// double leftBoundaryCondition(double y) {
//     return U_0 * (1 + alpha * V(y));
// }

// double rightBoundaryCondition(double y) {
//     return U_0;
// }

// double topBoundaryCondition(double x) {
//     return U_0;
// }

// double bottomBoundaryCondition(double x) {
//     return U_0;
// }

// // Fonction initiale devinée
// double initialGuessFunction(double x, double y) {
//     return U_0;
//     // return U_0 * (1 + alpha * V(y) * (1 - sin(M_PI * x / a)));
// }

// // Fonction du terme source

// double leftCoef = k_x * k_x;
// double rightCoef = k_y * k_y;
// double sourceTermFunction(double x, double y)
// {
//     // double leftTerm = ;
//     // double rightTerm = ;
//     return U_0 * alpha * (leftCoef * V(y) * sin(k_x * x) + rightCoef * V_tilde(x) * cos(k_y * y));
// }

// double exactSolutionFunction(double x, double y) {

//     return U_0 * (1. + alpha * V(y) * V_tilde(x));
//     // return leftBoundaryCondition(y);
// }

double a = 2.* M_PI;
double b = 2.* M_PI;
double U_0 = 0.;
// double alpha = 0.5;

double k_x = 2. * M_PI / a;
double k_y = 2. * M_PI / b;

double k_x_squared = k_x*k_x;
double k_y_squared = k_y*k_y;


double leftBoundaryCondition(double y) {
    return U_0;
}

double rightBoundaryCondition(double y) {
    return U_0;
}

double topBoundaryCondition(double x) {
    return U_0;
}

double bottomBoundaryCondition(double x) {
    return U_0;
}

double initialGuess(double x, double y) {
    return U_0;
}

double exactSolutionFunction(double x, double y) {
    return (double)sin(k_x * x) * sin(k_y * y);
}

// double leftCoef = (M_PI / (2 * a)) * (M_PI / (2 * a));
// double rightCoef = (2 * M_PI / b) * (2 * M_PI / b);
double sourceTermFunction(double x, double y)
{
    return (double)-(k_x_squared + k_y_squared) * exactSolutionFunction(x, y);
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    double x_min = 0.;
    double y_min = 0.;
    double x_max = a;
    double y_max = b;



    // int n = int(atoi(argv[1]));
    int n = 200;
    int Nx = n;
    int Ny = n;

    // printf("dx*dy = %.2g\n", d_x * d_y);
    // int maxIterations = 200000;
    int maxIterations = 20000;
    // int maxIterations = 20000;
    double epsilon = 1e-3;
    // std::cout << epsilon << "\n";

    // JacobiSequential()

    JacobiSequential jacobi_sequential(
        x_min, x_max, y_min, y_max, Nx, Ny,
        sourceTermFunction,
        leftBoundaryCondition,
        rightBoundaryCondition,
        topBoundaryCondition,
        bottomBoundaryCondition,
        initialGuess,
        exactSolutionFunction,
        true);
    // initialGuess
    jacobi_sequential.Solve(maxIterations, epsilon, true); // Max iterations: 1000, tolerance: 1e-6
    std::cout << "error = " << jacobi_sequential.error_ << "\n";

    // 363

    // JacobiMPI jacobi_mpi(
    //     x_min, x_max, y_min, y_max, Nx, Ny,
    //     sourceTermFunction,
    //     leftBoundaryCondition,
    //     rightBoundaryCondition,
    //     topBoundaryCondition,
    //     bottomBoundaryCondition,
    //     initialGuess,
    //     exactSolutionFunction,
    //     true
    //     );

    // jacobi_mpi.Solve(maxIterations, epsilon); // Max iterations: 1000, tolerance: 1e-6


    // std::cout << "error = " << jacobi_mpi.error_ << "\n";
    // 348
    MPI_Finalize();
    return 0;
}