MS01_Project_1
|— build/
|— visualization/
│	│— compare_fields.py
│	│— compare_logs.py
|— results/
|— utils/
│	│— Vector.hpp
│	│— Vector.hpp
|— solvers/
│	│— jacobi_sequential.hpp
│	│— jacobi_sequential.cpp
│	│— jacobi_MPI.cpp
│	│— jacobi_MPI.hpp
|— study_cases/
│	│— jacobi_sequential_asymptotic_accuracy/
│	│	│— main.cpp
│	│	│— results/
│	│— jacobi_MPI_asymptotic_accuracy/
│	│	│— main.cpp
│	│	│— results/
│	│— jacobi_sequential_MPI_comparison/
│	│	│— main.cpp
│	│	│— results/
|— main.cpp
|— CMakeLists.txt







double a = 20.;
double b = 2.;
double U_0 = 0.1;
double alpha = 0.5;

double k_y = 2. * M_PI / b;
double k_x = M_PI / (2. * a);
double V(double y)
{
    return 1 - cos(k_y * y);
}

double V_tilde(double x)
{
    return 1 - sin(k_x * x);
}

double leftBoundaryConditionFunction(double y) {
    return U_0 * (1 + alpha * V(y));
}

double rightBoundaryConditionFunction(double y) {
    return U_0;
}

double topBoundaryConditionFunction(double x) {
    return U_0;
}

double bottomBoundaryConditionFunction(double x) {
    return U_0;
}

// Fonction initiale devinée
double initialGuessFunction(double x, double y) {
    return U_0;
    // return U_0 * (1 + alpha * V(y) * (1 - sin(M_PI * x / a)));
}

// Fonction du terme source

double leftCoef = k_x * k_x;
double rightCoef = k_y * k_y;
double sourceTermFunction(double x, double y)
{
    // double leftTerm = ;
    // double rightTerm = ;
    return U_0 * alpha * (leftCoef * V(y) * sin(k_x * x) + rightCoef * V_tilde(x) * cos(k_y * y));
}

double exactSolutionFunction(double x, double y) {

    return U_0 * (1. + alpha * V(y) * V_tilde(x));
    // return leftBoundaryConditionFunction(y);
}
