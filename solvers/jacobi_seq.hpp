#ifndef JACOBI_SEQ_HPP
#define JACOBI_SEQ_HPP

#include "Vector.hpp"
#include "constants.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

class JacobiSeq {
public:
    int Nx_;
    int Ny_;
    int num_points_;
    double x_min_;
    double x_max_;
    double y_min_;
    double y_max_;
    double dx_;
    double dy_;

    int last_iteration_;    //! last iteration of the algorithm. (-1 if no convergence)

    Vector x_;       //! solution at iteration l
    Vector new_x_;   //! solution at iteration l + 1
    Vector r_;      //! residual vector at iteration l + 1 (r^(l + 1) = N(x^(l+1) - x^(l)) and r^(0) = b - Ax^(0))
    double r_0_;    //! initial residual (L2 norm of the initial residual vector)
    Vector x_exact_; //! exact solution if the exact solution function has been provided in the class constructor
    Vector b_;  //! source term/right-hand side vector

    // std::string algorithm_;     //! name of the algorithm (useful for the saved files)
    // std::string save_folder_;   //! folder path where to save the files (useful for the saved files)

    JacobiSeq( 
        const int Nx, const int Ny,
        const double x_min, const double x_max, 
        const double y_min, const double y_max,
        double (*source_term_function)(double, double),
        double (*left_boundary_condition_function)(double),
        double (*right_boundary_condition_function)(double),
        double (*top_boundary_condition_function)(double),
        double (*bottom_boundary_condition_function)(double),
        double (*initial_guess_function)(double, double),
        double (*exact_solution_function)(double, double) = nullptr);

    int I(int i, int j);

    void InitializeFields(
        double (*source_term_function)(double, double),
        double (*left_boundary_condition_function)(double),
        double (*right_boundary_condition_function)(double),
        double (*top_boundary_condition_function)(double),
        double (*bottom_boundary_condition_function)(double),
        double (*initial_guess_function)(double, double),
        double (*exact_solution_function)(double, double));

    void ComputeOneStep();

    void Solve(int max_iter, double epsilon);
};

// int Nx_;
// int Ny_;
// int num_points_;
// double x_min_;
// double x_max_;
// double y_min_;
// double y_max_;
// double dx_;
// double dy_;
// int last_iteration_;

#endif // JACOBI_SEQ_HPP