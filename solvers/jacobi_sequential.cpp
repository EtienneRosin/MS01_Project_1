#include "jacobi_sequential.hpp"

JacobiSequential::JacobiSequential(
    const double xMin, const double xMax, 
    const double yMin, const double yMax, 
    const int nX, const int nY,
    double (*sourceTermFunction)(double, double),
    double (*leftBoundaryConditionFunction)(double),
    double (*rightBoundaryConditionFunction)(double),
    double (*topBoundaryConditionFunction)(double),
    double (*bottomBoundaryConditionFunction)(double),
    double (*initialGuessFunction)(double, double),
    double (*exactSolutionFunction)(double, double),
    bool save_initialization, 
    std::string save_folder) : 
    Nx_(nX), Ny_(nY), num_points_((nX + 2)*(nY + 2)), 
    x_min_(xMin), x_max_(xMax), y_min_(yMin), y_max_(yMax), 
    // dx_((xMax - xMin)/(nX + 2.)), dy_((yMax - yMin)/(nY + 2.)), 
    dx_((xMax - xMin)/(nX + 1)), dy_((yMax - yMin)/(nY + 1)), 
    save_folder_(save_folder), last_iteration_(-1), 
    sol_(num_points_, 0.0), new_sol_(num_points_, 0.0), 
    source_(num_points_, 0.0), algorithm_(__func__){

    // std::cout << algorithm_ << "\n";
    // a_ii_ = -2. * ( (1./(dx_*dx_)) + (1./(dy_*dy_)) );
    // std::cout << "dx_ = " << dx_ << ", dy_ = "<< dy_ << "\n";

    if (exactSolutionFunction != nullptr)
    {
        exact_sol_.resize(num_points_, 0.0);
    }

    // Initializing the fields
    InitializeFields(
        sourceTermFunction,
        leftBoundaryConditionFunction,
        rightBoundaryConditionFunction,
        topBoundaryConditionFunction,
        bottomBoundaryConditionFunction,
        initialGuessFunction,
        exactSolutionFunction,
        save_initialization);

    if (exactSolutionFunction != nullptr)
    {
        ComputeError();
        std::cout << "initial error : " << error_ << "\n";
    }
    

    // ComputeResidual();
    initial_residual_ = std::sqrt(ComputeSquaredResidual(sol_, source_));
    // initial_residual_ = residual_;
    std::cout << "residu initi : " << initial_residual_ << "\n";
    
};

void JacobiSequential::InitializeFields(
    double (*sourceTermFunction)(double, double),
    double (*leftBoundaryConditionFunction)(double),
    double (*rightBoundaryConditionFunction)(double),
    double (*topBoundaryConditionFunction)(double),
    double (*bottomBoundaryConditionFunction)(double),
    double (*initialGuessFunction)(double, double),
    double (*exactSolutionFunction)(double, double),
    bool save_initialization) 
{
    // printf("Initializing fields...");
    for (int i = 0; i <= Nx_ + 1; i++) {
        double x = x_min_ + i * dx_;

        for (int j = 0; j <= Ny_ + 1; j++) {
            double y = y_min_ + j * dy_;

            source_[I(i, j)] = sourceTermFunction(x, y);
            if (exactSolutionFunction != nullptr)
            {
                exact_sol_[I(i, j)] = exactSolutionFunction(x, y);
            }

            if (i == 0) {
                // left boundary
                sol_[I(i, j)] = new_sol_[I(i, j)] = leftBoundaryConditionFunction(y);
            }
            else if (i == Nx_ + 1) {
                // right boundary
                // std::cout << "x_max : " << x_max_ << ", last x : " << x << "\n";
                sol_[I(i, j)] = new_sol_[I(i, j)] = rightBoundaryConditionFunction(y);
            }
            else if (j == 0) {
                // bottom boundary
                sol_[I(i, j)] = new_sol_[I(i, j)] = bottomBoundaryConditionFunction(x);
            }
            else if (j == Ny_ + 1) {
                // top boundary
                // std::cout << "y_max : " << y_max_ << ", last y : " << y << "\n";
                sol_[I(i, j)] = new_sol_[I(i, j)] = topBoundaryConditionFunction(x);
            }
            else {  
                // within the domain
                sol_[I(i, j)] = new_sol_[I(i, j)] = initialGuessFunction(x, y);
            }
        }
    }
    // std::cout << std::endl;
    if (save_initialization) {
        SaveOnDomain(sol_, "initial_guess");
        SaveOnDomain(source_, "source_term");
    if (exactSolutionFunction != nullptr)
    {
        SaveOnDomain(exact_sol_, "exact_solution");
    }
    }
    
    // std::cout.flush();
    // std::cout << "\nFields initialized.\n";
};


int JacobiSequential::I(int i, int j){
    // return i + (j - 1) * Nx_;
    return i + (Nx_ + 2) * j;
};

// void JacobiSequential::ComputeResidualVector()
double JacobiSequential::ComputeSquaredResidual(const Vector &u, const Vector &b) {
    double dx2 = dx_ * dx_;
    double dy2 = dy_ * dy_;
    // printf("On calcule le residu\n");
    double res_squared = 0.0;
    for (int i = 1; i <= Nx_; i++)
    {
        for (int j = 1; j <= Ny_; j++) {
            double tmp = b[I(i, j)];
            // std::cout << "b[I(i, j)] : " << b[I(i, j)] << "\n";
            // if (i == 1){
            //     std::cout << "value on bottom boundary : " << u[I(i - 1, j)] << "\n";
            // }
            // if (i == Nx_){
            //     std::cout << "value on top boundary : " << u[I(i + 1, j)] << "\n";
            // }
            // if (j == 1){
            //     std::cout << "value on left boundary : " << u[I(i, j-1)] << "\n";
            // }
            // if (j == Ny_){
            //     std::cout << "value on right boundary : " << u[I(i, j+1)] << "\n";
            // }
            tmp -= (u[I(i + 1, j)] - 2. * u[I(i, j)] + u[I(i - 1, j)]) / dx2;
            tmp -= (u[I(i, j + 1)] - 2. * u[I(i, j)] + u[I(i, j - 1)]) / dy2;

            res_squared += tmp * tmp;
        }
    }
    // printf("residu calculé\n");
    // std::cout << "res_squared = " << res_squared << "\n";
    return res_squared;
    // residual_ = std::sqrt(res_squared);

    // std::cout << "residual_ = " << residual_ << "\n";
};

void JacobiSequential::ComputeOneStep() {
    double dx2 = dx_ * dx_;
    double dy2 = dy_ * dy_;
    double A_ii = - 2. * (1. / dx2 + 1. / dy2);

    for (int i = 1; i <= Nx_; i++)
    {
        for (int j = 1; j <= Ny_; j++) {
            new_sol_[I(i, j)] = source_[I(i, j)] / A_ii;
            new_sol_[I(i, j)] -= (sol_[I(i + 1, j)] + sol_[I(i - 1, j)]) / (A_ii * dx2);
            new_sol_[I(i, j)] -= (sol_[I(i, j + 1)] + sol_[I(i, j - 1)]) / (A_ii * dy2);
        }
    }
    // std::cout << new_sol_ << "\n\n\n";
    // new_sol_.swap(sol_);
};


void JacobiSequential::Solve(int max_iter, double epsilon, bool save_solution) {
    std::ofstream log_file = CreateLogFile(max_iter, epsilon);
    // std::cout << "epsilon * initial_residual_ : " << epsilon * initial_residual_ << "\n";
    
    std::cout << "Computing..";
    for (int l = 1; l <= max_iter; l++) {
        ComputeOneStep();
        sol_.swap(new_sol_);
        // ComputeResidual();
        residual_ = (double)std::sqrt(ComputeSquaredResidual(sol_, source_));

        if (!exact_sol_.empty()) {
            error_ = (double)std::sqrt(ComputeSquaredError(sol_, exact_sol_));
            // ComputeError();
            log_file << l << "; " << residual_ << "; " << error_ << "\n";
            // std::cout << l << "; " << residual_ << "; " << error_ << "\n";
        }
        else {
            log_file << l << "; " << residual_ << "\n";
            // std::cout << l << "; " << residual_ << "\n";
        }
        
        // printf("Current iteration  \n");
        // std::cout << "l : " << l << ", residual_ : " << residual_ << "\n";

        // double diff = ComputeSquaredError(sol_, new_sol_);
        // if (std::sqrt(diff) <= epsilon){
        //     std::cout << "\nConverged after " << l << " iterations.\n";
        //     last_iteration_ = l;
        //     break;
        // }

        if (residual_ <= epsilon * initial_residual_) {
            std::cout.flush();
            std::cout << "\nConverged after " << l << " iterations.\n";
            last_iteration_ = l;
            break;
        }
    }
    
    std::cout.flush();
    std::cout << "Computing done.\n";
    log_file.close();

    if (save_solution){
        SaveOnDomain(sol_, "last_iteration_solution");
    }
    
};

void JacobiSequential::SaveOnDomain(const Vector &field, std::string file_name, std::string format) {
    std::string full_name = save_folder_ + algorithm_ + "_" + file_name + format;
    // printf("Saving : %s\n", full_name.c_str());
    std::ofstream file;
    file.open(full_name);
    file << "x; y; u\n";
    for (int i = 0; i <= Nx_ + 1; i++) {
        for (int j = 0; j <= Ny_ + 1; j++) {
            file << x_min_ + i * dx_ << "; " << y_min_ + j * dy_ << "; " << field[I(i, j)] << "\n";
        }
    }
    file.close();
};



double JacobiSequential::ComputeSquaredError(const Vector &u, const Vector &u_exact){
    double error_squared = 0.0;
    // ouble tmp = 0.0;
    for (int i = 0; i <= Nx_ + 1; i++)
    {
        for (int j = 0; j <= Ny_ + 1; j++){
            double tmp = exact_sol_[I(i, j)] - sol_[I(i, j)];
            error_squared += tmp * tmp;
        }
    }
    return error_squared;
}

void JacobiSequential::ComputeError() {
    error_ = 0.0;
    // Vector sub_solution = sol_.subVector(I(1, 1), I(Nx_ + 1, Ny_ + 1));
    // printf("exact_solution.size() = %zu, sub_solution.size() = %zu\n", exact_sol_.size(), sub_solution.size());
    // if (exact_sol_.size() != sub_solution.size())
    // {
    //     throw std::invalid_argument("Vector sizes do not match for subtraction");
    // }
    double tmp = 0.0;
    for (int i = 0; i <= Nx_ + 1; i++)
    {
        for (int j = 0; j <= Ny_ + 1; j++){
            double val = exact_sol_[I(i, j)] - sol_[I(i, j)];
            tmp += val * val;
        }
    }
    error_ = std::sqrt(tmp);
    // error_ = (exact_sol_ - sol_).norm();
}

// std::ofstream JacobiSequential::CreateLogFile(int max_iter, double epsilon, std::string format) {
    
//     std::string full_name = save_folder_ + algorithm_ + "_log" + format;
//     std::ofstream log_file;
//     log_file.open(full_name);

//     log_file << "dx; dy; l_max; epsilon\n";
//     log_file << dx_ << "; " << dy_<< "; " << max_iter << "; " << epsilon << "\n";

//     if (!exact_sol_.empty()) {
//         log_file << "Iteration; Residual; Error\n";
//     }
//     else {
//         log_file << "Iteration; Residual\n";
//     }
//     return log_file;
// };
std::ofstream JacobiSequential::CreateLogFile(int max_iter, double epsilon, std::string format) {
    std::string full_name = save_folder_ + algorithm_ + "_log" + format;
    std::ofstream log_file(full_name);
    
    if (!log_file.is_open()) {
        throw std::runtime_error("Could not open log file: " + full_name);
    }

    log_file << "dx; dy; l_max; epsilon\n";
    log_file << dx_ << "; " << dy_ << "; " << max_iter << "; " << epsilon << "\n";
    log_file << (exact_sol_.empty() ? "Iteration; Residual\n" : "Iteration; Residual; Error\n");
    
    return log_file;
}
// void Log(const std::string& message, std::ofstream& log_file) {
//     log_file << message << "\n"; // Can be extended for better formatting
// }