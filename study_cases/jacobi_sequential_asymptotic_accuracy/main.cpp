#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include "jacobi_sequential.hpp"

double a = 2. * M_PI;
double b = 2. * M_PI;
double U_0 = 1.;
// double alpha = 0.5;

double k_x = 2. * M_PI / a;
double k_y = 2. * M_PI / b;

double k_x_squared = pow(k_x, 2.);
double k_y_squared = pow(k_y, 2.);


double leftBoundaryCondition(double y) {
    return 0.0;
}

double rightBoundaryCondition(double y) {
    return 0.0;
}

double topBoundaryCondition(double x) {
    return 0.0;
}

double bottomBoundaryCondition(double x) {
    return 0.0;
}

double initialGuess(double x, double y) {
    return U_0;
}

double exactSolutionFunction(double x, double y) {
    return sin(k_x * x) * sin(k_y * y);
}

// double leftCoef = (M_PI / (2 * a)) * (M_PI / (2 * a));
// double rightCoef = (2 * M_PI / b) * (2 * M_PI / b);
double sourceTermFunction(double x, double y)
{
    return -(k_x_squared + k_y_squared) * exactSolutionFunction(x, y);
}



int main() {
    // Création du répertoire des résultats si nécessaire
    if (!std::filesystem::exists("results")) {
        std::filesystem::create_directory("results");
    }

    double xMin = 0.;
    double yMin = 0.;
    double xMax = a;
    double yMax = b;

    int n_min = 50; 
    int n_max = 1000; 
    double epsilon_min = 1e-1; 
    double epsilon_max = 1e-3;
    int maxIterations = 100000; 
    std::filesystem::path full_name = "../study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy.csv";
    // std::string full_name = std::filesystem::current_path().c_str() + "/results/jacobi_sequential_asymptotic_accuracy.csv";
    std::ofstream log_file(full_name);

    if (!log_file.is_open()) {
        
        std::cerr << "Error while opening the log file" << std::endl;
        return 1;
    }

    log_file << "epsilon; h; error; last_iteration\n";

    for (double epsilon = epsilon_min; epsilon >= epsilon_max; epsilon /= 10) {
        for (int n = n_min; n <= n_max; n += 100) {
            double h = (xMax - xMin) / (n + 1);
            std::cout << "epsilon = " << epsilon << ", h = " << h << "\n";
            int nX = n;
            int nY = n;

            JacobiSequential jacobi_sequential(
                xMin, xMax, yMin, yMax, nX, nY,
                sourceTermFunction,
                leftBoundaryCondition,
                rightBoundaryCondition,
                topBoundaryCondition,
                bottomBoundaryCondition,
                initialGuess,
                exactSolutionFunction,
                false,
                "../study_cases/jacobi_sequential_asymptotic_accuracy/results/"
            );

            jacobi_sequential.Solve(maxIterations, epsilon, false);
            log_file << epsilon << "; " << h << "; " << jacobi_sequential.error_ << "; " << jacobi_sequential.last_iteration_ << "\n";
        }
    }

    log_file.close();

    return 0;
}