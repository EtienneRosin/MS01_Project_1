#ifndef JACOBI_SEQUENTIAL_HPP
#define JACOBI_SEQUENTIAL_HPP

#include "Vector.hpp"
#include "constants.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>


class JacobiSequential {
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
    int last_iteration_;

    Vector sol_;
    Vector new_sol_;
    Vector source_;
    Vector exact_sol_;
    
    double initial_residual_;
    double residual_;
    
    double error_;
    
    // Vector residual_vector_;
    std::string algorithm_;
    std::string save_folder_;
    // double a_ii_;

    /*!
      \brief Default constructor is deleted to prevent its use.
    */
    JacobiSequential() = delete;

    /*!
      \brief Copy constructor is deleted to prevent copying of the JacobiMPI object.
    */
    JacobiSequential(const JacobiSequential &) = delete;

    /*!
      \brief Assignment operator is deleted to prevent assignment of the JacobiMPI object.
    */
    JacobiSequential &operator=(const JacobiSequential &) = delete;

    /*!
      \brief Constructor to initialize the JacobiSequential solver with the grid parameters, boundary conditions, and source term.

      \param xMin Minimum value of the x-axis.
      \param xMax Maximum value of the x-axis.
      \param yMin Minimum value of the y-axis.
      \param yMax Maximum value of the y-axis.
      \param nX Number of grid points along the x-axis.
      \param nY Number of grid points along the y-axis.
      \param sourceTermFunction Function that defines the source term.
      \param leftBoundaryConditionFunction Function that defines the left boundary condition.
      \param rightBoundaryConditionFunction Function that defines the right boundary condition.
      \param topBoundaryConditionFunction Function that defines the top boundary condition.
      \param bottomBoundaryConditionFunction Function that defines the bottom boundary condition.
      \param initialGuessFunction Function that provides the initial guess for the solution.
      \param exactSolutionFunction (Optional) Function that provides the exact solution for error calculation.
      \param save_initialization (Optional) If true, save the initial fields
      \param save_folder (Optional) folder where to save the files
    */
    JacobiSequential(
        const double xMin, const double xMax, 
        const double yMin, const double yMax, 
        const int nX, const int nY,
        double (*sourceTermFunction)(double, double),
        double (*leftBoundaryConditionFunction)(double),
        double (*rightBoundaryConditionFunction)(double),
        double (*topBoundaryConditionFunction)(double),
        double (*bottomBoundaryConditionFunction)(double),
        double (*initialGuessFunction)(double, double),
        double (*exactSolutionFunction)(double, double) = nullptr,
        bool save_initialization = false,
        std::string save_folder = "../results/");

    /*!
      \brief Initialize the following fields :
        - solution u
        - source term f
        - residual vector r
        - exact solution (if the exact solution function has been provided in the class constructor)
      \param sourceTermFunction Function that defines the source term.
      \param leftBoundaryConditionFunction Function that defines the left boundary condition.
      \param rightBoundaryConditionFunction Function that defines the right boundary condition.
      \param topBoundaryConditionFunction Function that defines the top boundary condition.
      \param bottomBoundaryConditionFunction Function that defines the bottom boundary condition.
      \param initialGuessFunction Function that provides the initial guess for the solution.
      \param exactSolutionFunction (Optional) Function that provides the exact solution for error calculation.
      \param save_initialization (Optional) If true, save the initial fields
    */
    void InitializeFields(
        double (*sourceTermFunction)(double, double),
        double (*leftBoundaryConditionFunction)(double),
        double (*rightBoundaryConditionFunction)(double),
        double (*topBoundaryConditionFunction)(double),
        double (*bottomBoundaryConditionFunction)(double),
        double (*initialGuessFunction)(double, double),
        double (*exactSolutionFunction)(double, double) = nullptr,
        bool save_initialization = false
    );

    /*!
      \brief Computes the linear index for a 2D grid based on the row and column indices.

      \param i Row index.
      \param j Column index.

      \return The linear index corresponding to the 2D grid point.
    */
    int I(int i, int j);

    /*!
      \brief Computes the error norm between the exact solution and the current approximation.
    */
    void ComputeError();
    double ComputeSquaredError(const Vector &u, const Vector &u_exact);

    /*!
      \brief Computes the local residual vector for the current solution.
      \todo Use the fact that the residual is given by : r^(l + 1) = N(x^(l + 1) - x^(l)). That may be quicker than the following implementation
    */
    void ComputeResidual();

    double ComputeSquaredResidual(const Vector &u, const Vector &b);

    /*!
      \brief Performs a single iteration of the Jacobi algorithm to update the solution.
    */
    void ComputeOneStep();

    /*!
      \brief Solves the system using the Jacobi alogrithm.

      \param max_iter Maximum number of iterations allowed (default is 1000).
      \param epsilon Tolerance for convergence (default is 1e-6).

      \details
      This function solves the system iteratively using the JacobiSequential method. At each iteration, the solution is updated and 
      the residual is computed to check if the method has converged. If the relative residual falls below the tolerance, 
      the iteration stops.
    */
    void Solve(int max_iter = kMaxIter, double epsilon = kEpsilon, bool save_solution = false);

    /*!
      \brief Save (in a txt file) the selected field values on the global domain
      \param field the considered field.
      \param file_name the file name.
      \param format file format (e.g. txt, csv...)
    */
    void SaveOnDomain(const Vector &field, std::string file_name, std::string format = ".csv");

    /*!
      \brief Create the log file containing the residuals and the errors (if an exact solution function has been given in the class constructor)
      \param max_iter Maximum number of iterations allowed (default is 1000).
      \param epsilon Tolerance for convergence (default is 1e-6).
      \param format file format (e.g. txt, csv...)

      \return the ofstream object associated to the file
      \warning you should close the stream after using it
    */
    std::ofstream CreateLogFile(int max_iter, double epsilon, std::string format = ".csv");

    
};

#endif // JACOBI_SEQUENTIAL_HPP