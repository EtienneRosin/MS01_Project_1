#include "jacobi_MPI.hpp"


JacobiMPI::JacobiMPI(
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
    bool save_initialization
    ) : global_Nx_(nX), global_Ny_(nY), x_min_(xMin), x_max_(xMax), y_min_(yMin), y_max_(yMax), 
    d_x_((xMax - xMin)/(nX + 1)), d_y_((yMax - yMin)/(nY + 1)) {

    algorithm_ = __func__;
    // Initialization of MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);

    // Getting rows indices of the subdomain
    int rows_per_process = nY / num_processes_ + 1;
    j_start_ = process_rank_ * rows_per_process + 1;                // First row of the subdomain
    j_stop_ = (process_rank_ + 1) * rows_per_process;               // Last row of the subdomain
    j_stop_ = (j_stop_ <= nY) ? j_stop_ : nY;                       // if j_stop > nY => j_stop = nY
    num_local_rows_ = j_stop_ - j_start_ + 1;                       // number of rows of the subdomain
    num_local_points_ = (global_Nx_ + 2) * (num_local_rows_ + 2);   // number of points in the subdomain (+2 for the ghost cells)
    // printf("p = %d, j_start_ =  %d, j_stop_ = %d\n", process_rank_, j_start_, j_stop_);
    // printf("num_local_rows_ = %d\n", num_local_rows_);

    // Resizing the fields
    local_solution_.resize(num_local_points_, 0.0);
    new_local_solution_.resize(num_local_points_, 0.0);
    local_source_term_.resize(num_local_points_, 0.0);
    local_residual_vector_.resize(num_local_points_, 0.0);

    if (process_rank_ == 0) {
        solution_.resize((global_Nx_ + 2) * (nY + 2));
    }

    if (exactSolutionFunction != nullptr)
    {
        local_exact_solution_.resize(num_local_points_, 0.0);
    }

    // Initializing the fields
    InitializeLocalFields(
        sourceTermFunction,
        leftBoundaryConditionFunction,
        rightBoundaryConditionFunction,
        topBoundaryConditionFunction,
        bottomBoundaryConditionFunction,
        initialGuessFunction,
        exactSolutionFunction,
        save_initialization);

    // Compute the initial residual
    ComputeLocalResidualVector();
    local_squared_residual_ = local_residual_vector_.NormSquared();

    CheckMPIError(
        MPI_Reduce(&local_squared_residual_, &global_squared_residual_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD),
        "Error reducing global squared residual");
    
    if (process_rank_ == 0) {
        initial_residual_ = std::sqrt(global_squared_residual_);
    }
        // initial_residual_ = std::sqrt(global_squared_residual_);
};


void JacobiMPI::InitializeLocalFields(
    double (*sourceTermFunction)(double, double),
    double (*leftBoundaryConditionFunction)(double),
    double (*rightBoundaryConditionFunction)(double),
    double (*topBoundaryConditionFunction)(double),
    double (*bottomBoundaryConditionFunction)(double),
    double (*initialGuessFunction)(double, double),
    double (*exactSolutionFunction)(double, double),
    bool save_initialization)
{
    for (int i = 0; i <= global_Nx_ + 1; i++) {
        double x = x_min_ + i * d_x_;
        if (process_rank_ == 0){    // top and bottom boundary condition of the global solution
            solution_[I(i, global_Ny_ + 1)] = topBoundaryConditionFunction(x);
            solution_[I(i, 0)] = bottomBoundaryConditionFunction(x);
        }

        for (int j = 0; j <= num_local_rows_ + 1; j++) {
            double y = y_min_ + (j + j_start_ - 1) * d_y_ ;
            
            if (i == 0) {
                // left boundary
                local_solution_[I(i, j)] = new_local_solution_[I(i, j)] = leftBoundaryConditionFunction(y);
            }
            else if (i == global_Nx_ + 1) {
                // right boundary
                local_solution_[I(i, j)] = new_local_solution_[I(i, j)] = rightBoundaryConditionFunction(y);
            }
            // ((j == 0) and (process_rank_ == 0))
            else if (j == 0) {
                // bottom boundary
                local_solution_[I(i, j)] = new_local_solution_[I(i, j)] = bottomBoundaryConditionFunction(x);
            }
            // ((j == num_local_rows_ + 1) and (process_rank_ == num_processes_ - 1))
            else if (j == num_local_rows_ + 1) {
                // top boundary
                local_solution_[I(i, j)] = new_local_solution_[I(i, j)] = topBoundaryConditionFunction(x);
            }
            else {  
                // within the domain
                local_solution_[I(i, j)] = initialGuessFunction(x, y);
                local_source_term_[I(i, j)] = sourceTermFunction(x, y);
                if (exactSolutionFunction != nullptr)
                {
                    local_exact_solution_[I(i, j)] = exactSolutionFunction(x, y);
                }
            }
        }
    }
    if (save_initialization) {
        SaveOnSubdomain(local_solution_, "initial_guess");
        SaveOnSubdomain(local_source_term_, "source_term");
        if (exactSolutionFunction != nullptr)
        {
            SaveOnSubdomain(local_exact_solution_, "exact_solution");
        }
    }
    
};

int JacobiMPI::I(int i, int j){
    return i + (global_Nx_ + 2) * j;
};

/*!
    \todo Use the fact that the residual is given by : r^(l + 1) = N(x^(l + 1) - x^(l)). That may be quicker than the following implementation
*/
void JacobiMPI::ComputeLocalResidualVector() {
    local_residual_vector_ = local_source_term_;
    for (int i = 1; i <= global_Nx_; i++) {
        for (int j = 1; j <= num_local_rows_; j++) {
            local_residual_vector_[I(i, j)] -= (local_solution_[I(i + 1, j)] - 2. * local_solution_[I(i, j)] + local_solution_[I(i - 1, j)]) / (d_x_ * d_x_);
            local_residual_vector_[I(i, j)] -= (local_solution_[I(i, j + 1)] - 2. * local_solution_[I(i, j)] + local_solution_[I(i, j - 1)]) / (d_y_ * d_y_);
        }
    }
};

void JacobiMPI::ComputeOneStep() {
    for (int i = 1; i <= global_Nx_; i++)
    {
        for (int j = 1; j <= num_local_rows_; j++) {
            double N_x_I = - 1. * (local_solution_[I(i + 1, j)] + local_solution_[I(i - 1, j)]) / (d_x_ * d_x_)
                           - 1. * (local_solution_[I(i, j + 1)] + local_solution_[I(i, j - 1)]) / (d_y_ * d_y_);
            double diag_coef = 1. / (d_x_ * d_x_) + 1. / (d_y_ * d_y_);
            new_local_solution_[I(i, j)] = -(1. / (2. * diag_coef)) * (N_x_I + local_source_term_[I(i, j)]);
        }
    }
};

void JacobiMPI::ExchangeBoundaries() {
    MPI_Request req_send_bottom, req_recv_top;
    if (process_rank_ > 0)
    { // => the subdomain has a neighbourd at its bottom
        // MPI_Isend(
        //     &local_solution_[I(1, 1)],  // first line
        //     global_Nx_,
        //     MPI_DOUBLE,
        //     process_rank_ - 1,
        //     0,
        //     MPI_COMM_WORLD,
        //     &req_send_bottom); // send the first line of the subdomain to its bottom neighbour
        // MPI_Irecv(
        //     &local_solution_[I(1, 0)],  // bottom ghost line
        //     global_Nx_,
        //     MPI_DOUBLE,
        //     process_rank_ - 1,
        //     0,
        //     MPI_COMM_WORLD,
        //     &req_recv_top); // receive the last line of the subdomain' bottom neighbour
        MPI_Isend(
            &local_solution_[I(0, 1)], // Include all points from left to right
            global_Nx_ + 2,            // Number of points along the x-direction
            MPI_DOUBLE,
            process_rank_ - 1,
            0,
            MPI_COMM_WORLD,
            &req_send_bottom);

        MPI_Irecv(
            &local_solution_[I(0, 0)], // Include ghost cells for the bottom row
            global_Nx_ + 2,
            MPI_DOUBLE,
            process_rank_ - 1,
            0,
            MPI_COMM_WORLD,
            &req_recv_top);
    }

    MPI_Request req_send_top, req_recv_bottom;
    if (process_rank_ < num_processes_ - 1)
    { // => the subdomain has a neighbourd at its top
        // MPI_Isend(
        //     // &local_solution_[I(1, num_local_rows_ + 1)],
        //     &local_solution_[I(1, num_local_rows_)], // top line
        //     global_Nx_,
        //     MPI_DOUBLE,
        //     process_rank_ + 1,
        //     0,
        //     MPI_COMM_WORLD,
        //     &req_send_top); // send the last line of the subdomain to its top neighbour
        // MPI_Irecv(
        //     // &local_solution_[I(1, num_local_rows_ + 2)],
        //     &local_solution_[I(1, num_local_rows_ + 1)],
        //     global_Nx_,
        //     MPI_DOUBLE,
        //     process_rank_ + 1,
        //     0,
        //     MPI_COMM_WORLD,
        //     &reqRecvTop); // receive the first line of the subdomain' top neighbour
        MPI_Isend(
            &local_solution_[I(0, num_local_rows_)], // top line
            global_Nx_ + 2,
            MPI_DOUBLE,
            process_rank_ + 1,
            0,
            MPI_COMM_WORLD,
            &req_send_top); // send the last line of the subdomain to its top neighbour
        MPI_Irecv(
            &local_solution_[I(0, num_local_rows_ + 1)],
            global_Nx_ + 2,
            MPI_DOUBLE,
            process_rank_ + 1,
            0,
            MPI_COMM_WORLD,
            &req_recv_bottom); // receive the first line of the subdomain' top neighbour
    }

    if (process_rank_ > 0) {
        
        MPI_Wait(&req_send_bottom, MPI_STATUS_IGNORE);
        MPI_Wait(&req_recv_top, MPI_STATUS_IGNORE);
    }
    if (process_rank_ < num_processes_ - 1) {
        MPI_Wait(&req_send_top, MPI_STATUS_IGNORE);
        MPI_Wait(&req_recv_bottom, MPI_STATUS_IGNORE);
    }
};

void JacobiMPI::Solve(int max_iter, double epsilon, bool save_solution) {
    bool converged = false;

    std::ofstream log_file;
    if (process_rank_ == 0) {
        log_file = CreateLogFile(max_iter, epsilon);
    }
    ExchangeBoundaries();
    for (int l = 0; l < max_iter; l++) {
        if (converged) break;
        // ExchangeBoundaries();
        ComputeOneStep();
        local_solution_.swap(new_local_solution_);
        ExchangeBoundaries();
        if (!local_exact_solution_.empty()) {
            ComputeErrorSquared();
            CheckMPIError(
            MPI_Reduce(&local_squared_error_, &global_squared_error_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD),
            "Error reducing global squared error");
        }

        // Check the residual stopping criteria
        ComputeLocalResidualVector();
        local_squared_residual_ = local_residual_vector_.NormSquared();

        CheckMPIError(
            MPI_Reduce(&local_squared_residual_, &global_squared_residual_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD),
            "Error reducing global squared residual");

        if (process_rank_ == 0) {
            double relative_residual = std::sqrt(global_squared_residual_) / initial_residual_;
            if (!local_exact_solution_.empty()) {
                log_file << l << "; " << relative_residual << "; " << std::sqrt(global_squared_error_) << "\n";
            }
            else {
                log_file << l << "; " << relative_residual << "\n";
            }

            if (relative_residual <= epsilon) {
                std::cout << "Converged after " << l << " iterations.\n";
                converged = true;
            }
        }
        CheckMPIError(
            MPI_Bcast(&converged, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD),
            "Error broadcasting convergence signal");
    }
    if (process_rank_ == 0) {
        log_file.close();
    }
    // Save the solution
    MPI_Barrier(MPI_COMM_WORLD);
    if (save_solution) {
        SaveOnSubdomain(local_solution_, "last_iteration_solution");
    }
    
    GatherGlobalSolution();
    if ((process_rank_ == 0) and save_solution) {
        SaveOnGlobalDomain(solution_, "last_iteration_solution");
    }
};

void JacobiMPI::GatherGlobalSolution() {
    MPI_Gather(
        &local_solution_[I(0, 1)],          // starting address of send buffer
        (global_Nx_ + 2) * num_local_rows_, // number of elements in send buffer
        MPI_DOUBLE,                         // data type of send buffer elements
        &solution_[I(0, 1)],         // Starting address of receive buffer
        (global_Nx_ + 2) * num_local_rows_, // number of elements for any single receive
        MPI_DOUBLE,                         // data type of recv buffer elements
        0,                                  // rank of receiving process
        MPI_COMM_WORLD);                    // communicator
};

void JacobiMPI::SaveOnSubdomain(const Vector &field, std::string file_name, std::string folder_path, std::string format) {
    std::string full_name = folder_path + algorithm_ + "_" + file_name + "_" + std::to_string(process_rank_) + format;
    // printf("Saving : %s\n", full_name.c_str());
    std::ofstream file;
    file.open(full_name);
    file << "x; y; u\n";
    for (int i = 0; i <= global_Nx_ + 1; i++) {
        for (int j = 0; j <= num_local_rows_ + 1; j++) {
            file << x_min_ + i * d_x_ << "; " << y_min_ + (j + j_start_ - 1) * d_y_ << "; " << field[I(i, j)] << "\n";
        }
    }
    file.close();
};

void JacobiMPI::SaveOnGlobalDomain(const Vector &field, std::string file_name, std::string folder_path, std::string format) {
    std::string full_name = folder_path + algorithm_ + "_" + file_name + format;
    std::ofstream file;
    file.open(full_name);
    file << "x; y; u\n";
    for (int i = 0; i <= global_Nx_ + 1; i++) {
        for (int j = 0; j <= global_Ny_ + 1; j++) {
            file << x_min_ + i * d_x_ << "; " << y_min_ + j * d_y_ << "; " << field[I(i, j)] << "\n";
        }
    }
    file.close();
};

void JacobiMPI::ComputeErrorSquared() {
    local_squared_error_ = 0.0;
    int j_min = 1;
    int j_max = num_local_rows_;
    
    // if (process_rank_ == 0) {
    //     j_min -= 1;
    // }
    // else if (process_rank_ == num_processes_ - 1){
    //     j_max += 1;
    // }

    for (int i = 1; i <= global_Nx_; i++)
    {
        for (int j = 1; j <= num_local_rows_; j++) {
            local_squared_error_ += pow(local_exact_solution_[I(i, j)] - local_solution_[I(i, j)], 2);
        }
    }
        // local_squared_error_ = (local_exact_solution_ - local_solution_).dot(local_exact_solution_ - local_solution_);
};

std::ofstream JacobiMPI::CreateLogFile(int max_iter, double epsilon, std::string folder_path, std::string format) {
    std::string full_name = folder_path + algorithm_ + "_log" + format;
    std::ofstream log_file;
    log_file.open(full_name);

    log_file << "dx; dy; l_max; epsilon\n";
    log_file << d_x_ << "; " << d_y_<< "; " << max_iter << "; " << epsilon << "\n";

    if (!local_exact_solution_.empty()) {
        log_file << "Iteration; Residual; Error\n";
    }
    else {
        log_file << "Iteration; Residual\n";
    }
    return log_file;
};


void JacobiMPI::CheckMPIError(int errCode, const std::string& errMsg) {
    if (errCode != MPI_SUCCESS) {
        std::cerr << "MPI Error on process " << process_rank_ << " : " << errMsg << std::endl;
        MPI_Abort(MPI_COMM_WORLD, errCode);
    }
};