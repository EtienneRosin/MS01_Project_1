#include "jacobi_seq.hpp"


JacobiSeq::JacobiSeq(
    const int Nx, const int Ny,
    const double x_min, const double x_max, 
    const double y_min, const double y_max,
    double (*source_term_function)(double, double),
    double (*left_boundary_condition_function)(double),
    double (*right_boundary_condition_function)(double),
    double (*top_boundary_condition_function)(double),
    double (*bottom_boundary_condition_function)(double),
    double (*initial_guess_function)(double, double),
    double (*exact_solution_function)(double, double)) : 
    Nx_(Nx), Ny_(Ny), num_points_((Nx + 2)*(Ny + 2)), 
    x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), 
    dx_((x_max - x_min)/(Nx + 1)), dy_((y_max - y_min)/(Ny + 1)), 
    last_iteration_(-1) {
    

    x_.resize(num_points_, 0.0);
    new_x_.resize(num_points_, 0.0);
    r_.resize(num_points_, 0.0);
    b_.resize(num_points_, 0.0);

    if (exact_solution_function != nullptr){
        x_exact_.resize(num_points_, 0.0);
    }

    InitializeFields(
        source_term_function,
        left_boundary_condition_function,
        right_boundary_condition_function,
        top_boundary_condition_function,
        bottom_boundary_condition_function,
        initial_guess_function,
        exact_solution_function);

    // compute initial residual


    // std::cout << r_ << "\n";
    for (int i = 1; i <= Nx_; i++) {
        for (int j = 1; j <= Ny_; j++) {
            r_[I(i, j)] = b_[I(i,j)] - ( ( x_[I(i + 1,j)] - 2. * x_[I(i,j)] + x_[I(i - 1,j)] )/pow(dx_, 2.) + (x_[I(i,j + 1)] - 2. * x_[I(i,j)] + x_[I(i,j - 1)] )/pow(dy_, 2.));
        }
    }
    // std::cout << r_ << "\n";
    r_0_ = r_.norm();
    // std::cout << "r_O_ = " << r_0_ << "\n";

    // Vector r_;      //! residual vector at iteration l + 1 (r^(l + 1) = N(x^(l+1) - x^(l)) and r^(0) = b - Ax^(0))
    // double r_0_;    //! initial residual (L2 norm of the initial residual vector)
    // Vector x_exact_; //! exact solution if the exact solution function has been provided in the class constructor
    // Vector b_;  //! source term/right-hand side vector
};

int JacobiSeq::I(int i, int j){
    // return i + (j - 1) * Nx_;
    return i + (Nx_ + 2) * j;
};


void JacobiSeq::InitializeFields(
    double (*sourceTermFunction)(double, double),
    double (*leftBoundaryConditionFunction)(double),
    double (*rightBoundaryConditionFunction)(double),
    double (*topBoundaryConditionFunction)(double),
    double (*bottomBoundaryConditionFunction)(double),
    double (*initialGuessFunction)(double, double),
    double (*exactSolutionFunction)(double, double)) 
{
    printf("Initializing fields...");
    for (int i = 0; i <= Nx_ + 1; i++) {
        double x = x_min_ + i * dx_;
        for (int j = 0; j <= Ny_ + 1; j++) {
            double y = y_min_ + j * dy_;

            if (i == 0) {
                // left boundary
                x_[I(i, j)] = new_x_[I(i, j)] = leftBoundaryConditionFunction(y);
            }
            else if (i == Nx_ + 1) {
                // right boundary
                x_[I(i, j)] = new_x_[I(i, j)] = rightBoundaryConditionFunction(y);
            }
            else if (j == 0) {
                // bottom boundary
                x_[I(i, j)] = new_x_[I(i, j)] = bottomBoundaryConditionFunction(x);
            }
            else if (j == Ny_ + 1) {
                // top boundary
                x_[I(i, j)] = new_x_[I(i, j)] = topBoundaryConditionFunction(x);
            }
            else {  
                // within the domain
                x_[I(i, j)] = new_x_[I(i, j)] = initialGuessFunction(x, y);
                // b_[I(i, j)] = sourceTermFunction(x, y);
                b_[I(i, j)] = sourceTermFunction(x, y);
            }
            
            if (exactSolutionFunction != nullptr)
            {
                x_exact_[I(i, j)] = exactSolutionFunction(x, y);
            }

            
        }
    }
    // std::cout << std::endl;
    // if (save_initialization) {
    //     SaveOnDomain(x_, "initial_guess");
    //     SaveOnDomain(b_, "source_term");
    // if (exactSolutionFunction != nullptr)
    // {
    //     SaveOnDomain(x_exact_, "exact_solution");
    // }
    // }
    
    std::cout.flush();
    std::cout << "\nFields initialized.\n";
};


void JacobiSeq::ComputeOneStep(){
    for (int i = 1; i <= Nx_; i++) {
        for (int j = 0; j <= Ny_; j++) {
            new_x_[I(i, j)] = b_[I(i, j)] - (x_[I(i + 1, j)] + x_[I(i - 1, j)]) / pow(dx_, 2.) - (x_[I(i, j + 1)] + x_[I(i, j - 1)]) / pow(dy_, 2.);

            r_[I(i, j)] = new_x_[I(i, j)] + 2 * (1 / pow(dx_, 2.) + 1 / pow(dy_, 2.)) * x_[I(i, j)];

            new_x_[I(i, j)] *= -1 / (2 * (1 / pow(dx_, 2.) + 1 / pow(dy_, 2.)));
            // r_[I(i, j)] = b_[I(i, j)] - ((x_[I(i + 1, j)] - 2. * x_[I(i, j)] + x_[I(i - 1, j)]) / pow(dx_, 2.) + (x_[I(i, j + 1)] - 2. * x_[I(i, j)] + x_[I(i, j - 1)]) / pow(dy_, 2.));
        }
    }
    x_.swap(new_x_);

    
};

void JacobiSeq::Solve(int max_iter, double epsilon) {
    std::cout << "Computing..";
    for (int l = 0; l < max_iter; l++) {
        ComputeOneStep();
        double relative_residual = r_.norm() / r_0_;
        // std::cout << "relative_residual = " << relative_residual << "\n";

        if (relative_residual <= epsilon) {
            std::cout.flush();
            std::cout << "\nConverged after " << l << " iterations.\n";

            if (!x_exact_.empty()){
                std::cout << "last iteration error : " << (x_exact_ - x_).norm() << "\n";
            }

            last_iteration_ = l;
            break;
        }
    }
    // std::cout << "\n No convergence\n";
    // if (!x_exact_.empty()){
    //     std::cout << "last iteration error : " << (x_exact_ - x_).norm() << "\n";
    // }
};