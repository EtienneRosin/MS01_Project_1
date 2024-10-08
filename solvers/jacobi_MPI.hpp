#ifndef JACOBI_MPI_HPP
#define JACOBI_MPI_HPP

#include "Vector.hpp"
#include "constants.hpp"

#include <mpi.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


class JacobiMPI {
public:
    int global_Nx_;
    int global_Ny_;
    int num_local_rows_;
    int num_local_points_;
    int j_start_;
    int j_stop_;
    double initial_residual_;
    double local_squared_residual_;
    double global_squared_residual_;
    double x_min_;
    double x_max_;
    double y_min_;
    double y_max_;
    double d_x_;
    double d_y_;
    int process_rank_;
    int num_processes_;
    double local_squared_error_;
    double global_squared_error_;
    Vector local_solution_;
    Vector new_local_solution_;
    Vector local_source_term_;
    Vector local_residual_vector_;
    Vector local_exact_solution_;
    Vector solution_;
    bool save_initialization;

    double a_ii_;
    double error_;

    std::string algorithm_;

    /*!
      \brief Default onstructor is deleted to prevent its use.
    */
    JacobiMPI() = delete;

    /*!
      \brief Copy constructor is deleted to prevent copying of the JacobiMPI object.
    */
    JacobiMPI(const JacobiMPI &) = delete;

    /*!
      \brief Assignment operator is deleted to prevent assignment of the JacobiMPI object.
    */
    JacobiMPI &operator=(const JacobiMPI &) = delete;

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
    */
    JacobiMPI(
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
        bool save_initialization = false);

    /*!
      \brief Initialize the following fields on the process' subdomain :
        - local solution u
        - local source term f
        - local residual vector r
        - local exact solution (if the exact solution function has been provided in the class constructor)
      \param sourceTermFunction Function that defines the source term.
      \param leftBoundaryConditionFunction Function that defines the left boundary condition.
      \param rightBoundaryConditionFunction Function that defines the right boundary condition.
      \param topBoundaryConditionFunction Function that defines the top boundary condition.
      \param bottomBoundaryConditionFunction Function that defines the bottom boundary condition.
      \param initialGuessFunction Function that provides the initial guess for the solution.
      \param exactSolutionFunction (Optional) Function that provides the exact solution for error calculation.
      \param save_initialization (Optional) If true, save the initial fields
    */
    void InitializeLocalFields(
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
      \brief Computes the error norm squared between the exact solution and the current approximation.
    */
    void ComputeErrorSquared();

    /*!
      \brief Computes the local residual vector for the current solution.
      \todo Use the fact that the residual is given by : r^(l + 1) = N(x^(l + 1) - x^(l)). That may be quicker than the following implementation
    */
    void ComputeLocalResidualVector();

    /*!
      \brief Performs a single iteration of the Jacobi algorithm to update the solution.
    */
    void ComputeOneStep();

    /*!
      \brief Solves the system using the Jacobi alogrithm.

      \param max_iter Maximum number of iterations allowed (default is 1000).
      \param epsilon Tolerance for convergence (default is 1e-6).
      \param save_solution (Optional) if true, save the last iteration solution

      \details
      This function solves the system iteratively using the JacobiSequential method. At each iteration, the solution is updated and 
      the residual is computed to check if the method has converged. If the relative residual falls below the tolerance, 
      the iteration stops.
    */
    void Solve(int max_iter = kMaxIter, double epsilon = kEpsilon, bool save_solution = false);

    /*!
      \brief Exchange the boundary values to the neighbour subdomain.
    */
    void ExchangeBoundaries();

    /*!
      \brief Gather all the local solutions
    */
    void GatherGlobalSolution();

    /*!
      \brief Save (in a txt file) the selected field values on the subdomain.
      \param field the considered field.
      \param file_name the file name (that would be completed with the process number).
      \param folder_path path of the saving folder.
      \param format file format (e.g. txt, csv...)
    */
    void SaveOnSubdomain(const Vector &field, std::string file_name, std::string folder_path = "../results/", std::string format = ".csv");

    /*!
      \brief Save (in a txt file) the selected field values on the global domain
      \param field the considered field.
      \param file_name the file name.
      \param folder_path path of the saving folder.
      \param format file format (e.g. txt, csv...)
    */
    void SaveOnGlobalDomain(const Vector &field, std::string file_name, std::string folder_path = "../results/", std::string format = ".csv");

    /*!
      \brief Create the log file containing the residuals and the errors (if an exact solution function has been given in the class constructor)
      \param max_iter Maximum number of iterations allowed (default is 1000).
      \param epsilon Tolerance for convergence (default is 1e-6).
      \param folder_path path of the saving folder.
      \param format file format (e.g. txt, csv...)

      \return the ofstream object associated to the file
      \warning you should close the stream after using it
    */
    std::ofstream CreateLogFile(int max_iter, double epsilon, std::string folder_path = "../results/", std::string format = ".csv");

    /*!
      \brief Check if the MPI function throws an error
      \param errCode error code
      \param errMsg error message
      \throw std::cerr : the error message
    */
    void CheckMPIError(int errCode, const std::string &errMsg);
};

#endif // JACOBI_MPI_HPP