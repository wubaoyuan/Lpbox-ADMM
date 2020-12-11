//
// Created by squintingsmile on 2020/9/22.
//

#include <cstring>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include "solver.h"
Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
                             ", ", ", ", "", "", " << ", ";");

int solver::conjugate_gradient_with_guess(const SparseMatrix &A, const DenseVector &b, const DenseVector &init_guess,
                                          DenseVector &result, double tol, int max_iter, bool printMsg) {
    DenseVector r = b - A * init_guess;
    if (r.norm() <= tol * b.norm()) {
        result = init_guess;
        printf("The initial guess is within the convergence tolerance");
        return 1;
    }
    DenseVector p = r;
    DenseVector Ap;
    DenseVector x = init_guess;
    double r0tr0 = 0;
    double r1tr1 = r.squaredNorm();
    double alpha = 0;
    double beta = 0;
    for (int i = 0; i < max_iter; i++) {
        Ap = A * p;
        alpha = r1tr1 / (p.transpose() * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        if (r.norm() <= tol * b.norm()) {
            result = x;
            if (printMsg)
                printf("Conjugate gradient converges after %d iterations\n", i);
            return 1;
        }
        r0tr0 = r1tr1;
        r1tr1 = r.squaredNorm();
        beta = r1tr1 / r0tr0;
        p = r + beta * p;
    }

    result = x;
    printf("Conjugate gradient does not converge after iteration number of %d\n", max_iter);
    return 1;
}

int solver::conjugate_gradient_with_guess_pre(const SparseMatrix &A, const DenseVector &b, const DenseVector &init_guess,
                                                      DenseVector &result, double tol, int max_iter,
                                                      const Preconditioner preconditioner, bool printMsg) {
    if (preconditioner == INCOMPLETE_CHOLESKY) {

    }
    return 1;
}


float_t solver::std_dev(std::vector<float_t> arr, size_t begin, size_t end) {
    float_t mean = 0;
    float_t std_deviation = 0;
    size_t size = end - begin;
    for (int i = begin; i < end; i++) {
        mean += arr[i];
    }

    mean /= size;
    for (int i = 0; i < size; i++) {
        std_deviation += (arr[begin + i] - mean) * (arr[begin + i] - mean);
    }

    std_deviation /= size - 1;
    if (std_deviation == 0) {
        return 0;
    }

    return std::pow(std_deviation, 1.0 / 2);
}

void solver::project_box(int n, const double *x, double *y) {
    y = new double[n];
    while (n--) {
        if (x[n] > 1) {
            y[n] = 1;
        } else {
            if (x[n] < 0) {
                y[n] = 0;
            } else {
                y[n] = x[n];
            }
        }
    }
}

void solver::project_box(int n, const DenseVector &x, DenseVector &y) {
    while (n--) {
        if (x[n] > 1) {
            y[n] = 1;
        } else {
            if (x[n] < 0) {
                y[n] = 0;
            } else {
                y[n] = x[n];
            }
        }
    }
}

void solver::project_shifted_Lp_ball(int n, const DenseVector &x, int p, DenseVector &y) {
    y.array() = x.array() - 0.5;
    float_t normp_shift = y.norm();

    y.array() = y.array() * std::pow(n, 1.0 / p) / (2 * normp_shift) + 0.5;
}

float_t solver::compute_cost(const DenseVector& x, const SparseMatrix& A, const DenseVector& b) {
    auto val = (x.transpose() * A).dot(x);
    auto val2 = b.dot(x);
    return val + val2;
}

float_t solver::compute_std_obj(std::vector<float_t> obj_list, int history_size) {
    size_t s = obj_list.size();
    float_t std_obj;
    if (s <= history_size) {
        std_obj = std_dev(obj_list, 0, s);
    } else {
        std_obj = std_dev(obj_list, s - history_size, s);
    }

    return std_obj / std::abs(obj_list[s-1]);
}

void solver::ADMM_bqp_unconstrained_init() {
    std_threshold = 1e-6;
    gamma_val = 1.0;
    gamma_factor = 0.99;
    initial_rho = 5;
    learning_fact = 1 + 3.0 / 100;
    rho_upper_limit = 1000;
    history_size = 5;
    rho_change_step = 5;
    rel_tol = 1e-5;
    stop_threshold = 1e-3;
    max_iters = 1e4;
    projection_lp = 2;
    pcg_tol = 1e-3;
    pcg_maxiters = 1e3;
}

void solver::ADMM_bqp_linear_eq_init() {
    stop_threshold = 1e-4;
    std_threshold = 1e-6;
    gamma_val = 1.6;
    gamma_factor = 0.95;
    rho_change_step = 5;
    max_iters = 1e3;
    initial_rho = 25;
    history_size = 3;
    learning_fact = 1 + 1.0 / 100;
    pcg_tol = 1e-4;
    pcg_maxiters = 1e3;
    rel_tol = 5e-5;
    projection_lp=2;
}

void solver::ADMM_bqp_linear_ineq_init() {
    stop_threshold = 1e-4;
    std_threshold = 1e-6;
    gamma_val = 1.6;
    gamma_factor = 0.95;
    rho_change_step = 5;
    max_iters = 1e3;
    initial_rho = 25;
    history_size = 3;
    learning_fact = 1 + 1.0 / 100;
    pcg_tol = 1e-4;
    pcg_maxiters = 1e3;
    rel_tol = 5e-5;
    projection_lp = 2;

}

void solver::ADMM_bqp_linear_eq_and_uneq_init() {
    stop_threshold = 1e-4;
    std_threshold = 1e-6;
    gamma_val = 1.6;
    gamma_factor = 0.95;
    rho_change_step = 5;
    max_iters = 1e3;
    initial_rho = 25;
    history_size = 3;
    learning_fact = 1 + 1.0 / 100;
    pcg_tol = 1e-4;
    pcg_maxiters = 1e3;
    rel_tol = 5e-5;
    projection_lp = 2;
}

/** This function solves the following optimization problem min_x (x^T A x + b^T x) such that x is {0,1}^n
 *  ADMM update steps with x0 is feasible and binary. Currently this function is using the double precision
 *  version of the BLAS and LAPACK library. The input matrix is column major matrix.
 * @param A: The matrix A
 * @param b: The vector b
 * @param n: The number of columns of A (also the size of vector b and x)
 */

int solver::ADMM_bqp_unconstrained(int matrix_layout, int n, const SparseMatrix &_A,
                                   const DenseVector &_b, Solution& sol, bool printMsg) {
    auto x_sol    = DenseVector(n);
    auto y1       = DenseVector(n);
    auto y2       = DenseVector(n);
    auto z1       = DenseVector(n);
    auto z2       = DenseVector(n);
    auto prev_idx = DenseVector(n);
    auto best_sol = DenseVector(n);
    auto temp_vec = DenseVector(n);
    auto _2A      = SparseMatrix(n, n);
    auto temp_mat = SparseMatrix(n, n);
    auto cur_idx  = DenseVector(n);
    float_t cur_obj;

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> cg;
    cg.setMaxIterations(pcg_maxiters);
    cg.setTolerance(pcg_tol);

    _2A = 2 * _A;

    if (x0 != nullptr) {
        memcpy(x_sol.data(), x0, sizeof(float_t) * n);
    } else {
        memset(x_sol.data(), 0, sizeof(float_t) * n);
    }

    float_t rho1 = initial_rho;
    float_t rho2 = initial_rho;
    std::vector<float_t> obj_list; /* Stores the objective value calculated during each iteration */
    float_t std_obj = 1;

    y1 = x_sol;
    y2 = x_sol;
    prev_idx = x_sol;
    best_sol = x_sol;

    float_t best_bin_obj = compute_cost(x_sol, _A, _b);

//    std::cout << _A << std::endl;

    if (printMsg) {
        std::cout << "Initial state" << std::endl;
        std::cout << "norm of x_sol: " << x_sol.norm() << std::endl;
        std::cout << "norm of y1: " << y1.norm() << std::endl;
        std::cout << "norm of y2: " << y2.norm() << std::endl;
        std::cout << "norm of z1: " << z1.norm() << std::endl;
        std::cout << "norm of z2: " << z2.norm() << std::endl;
        std::cout << "norm of cur_idx: " << cur_idx.norm() << std::endl;
        std::cout << "-------------------------------------------------" << std::endl;
    }

    long time_elapsed = 0;
    for (int iter = 0; iter < max_iters; iter++) {
//        std::cout << x_sol << std::endl;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        temp_vec = x_sol + z1 / rho1;

        /* Project vector on [0, 1] box */
        project_box(n, temp_vec, y1);

        temp_vec = x_sol + z2 / rho2;

        /* Project vector on shifted lp box */
        project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

        temp_mat = _2A;
        for (int i = 0; i < n; i++) {
            temp_mat.coeffRef(i, i) += rho1 + rho2;
        }
        temp_vec = rho1 * y1 + rho2 * y2 - (_b + z1 + z2);
        temp_mat.makeCompressed();
        /* The solveWithGuess function does not seems to check the tolerance condition in the beginning */
        if ((temp_mat * y1 - temp_vec).norm() < cg.tolerance() * temp_vec.norm()) {
            x_sol = y1;
        } else {
            cg.compute(temp_mat);
            x_sol = cg.solveWithGuess(temp_vec, y1);
//            conjugate_gradient_with_guess(temp_mat, temp_vec, y1, x_sol, pcg_tol, pcg_maxiters, true);
            auto info = cg.info();
            if (info != Eigen::Success) {
                std::cout << "conjugate gradient error info: " << info << std::endl;
                exit(-1);
            }
        }



        z1 = z1 + gamma_val * rho1 * (x_sol - y1);
        z2 = z2 + gamma_val * rho2 * (x_sol - y2);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time_elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

        if ((iter+2) % rho_change_step == 0) {
            rho1 = learning_fact * rho1;
            rho2 = learning_fact * rho2;
            gamma_val = std::max(gamma_val * gamma_factor, float_t(1.0));
            std::cout << "rho1: " <<  rho1 << std::endl;
        }

        if (iter == 10000) {
            //    std::ofstream temp_mat_c;
            //    std::ofstream temp_vec_c;
            std::ofstream temp_sol_c;
            std::ofstream temp_init_c;
            //    temp_mat_c.open("temp_mat_c");
            //    temp_vec_c.open("temp_vec_c");
            temp_sol_c.open(("temp_sol_c" + std::to_string(100)).c_str());
            temp_init_c.open(("temp_init_c" + std::to_string(100)).c_str());

            std::ofstream temp_y1_c;
            std::ofstream temp_y2_c;
            temp_y1_c.open(("temp_y1_c" + std::to_string(100)).c_str());
            temp_y2_c.open(("temp_y2_c" + std::to_string(100)).c_str());
            std::ofstream temp_z1_c;
            std::ofstream temp_z2_c;
            temp_z1_c.open(("temp_z1_c" + std::to_string(100)).c_str());
            temp_z2_c.open(("temp_z2_c" + std::to_string(100)).c_str());

            std::ofstream temp_rho1_c;
            std::ofstream temp_rho2_c;
            temp_rho1_c.open(("temp_rho1_c" + std::to_string(100)).c_str());
            temp_rho2_c.open(("temp_rho2_c" + std::to_string(100)).c_str());

            auto output_vec = Eigen::MatrixXd(temp_mat);
            output_vec.resize(n * n, 1);
//                temp_mat_c << output_vec;
//                temp_vec_c << temp_vec;
            temp_sol_c << x_sol;
            temp_init_c << y1;
            temp_init_c.close();
//                temp_mat_c.close();
//                temp_vec_c.close();
            temp_sol_c.close();
            temp_y1_c << y1;
            temp_y1_c.close();
            temp_y2_c << y2;
            temp_y2_c.close();
            temp_z1_c << z1;
            temp_z1_c.close();
            temp_z2_c << z2;
            temp_z2_c.close();
            temp_rho1_c << rho1;
            temp_rho1_c.close();
            temp_rho2_c << rho2;
            temp_rho2_c.close();

            exit(0);
        }

        float_t temp0 = std::max(x_sol.norm(), float_t(2.2204e-16));
        float_t temp1 = (x_sol-y1).norm() / temp0;
        float_t temp2 = (x_sol-y2).norm() / temp0;
        if (temp1 <= stop_threshold && temp2 <= stop_threshold) {
            printf("iter: %d, stop_threshold: %.6f\n", iter, std::max(temp1,temp2));
            break;
        }

        float_t obj_val = compute_cost(x_sol,_A,_b);
        obj_list.push_back(obj_val);
        if (obj_list.size() >= history_size) {
            std_obj = compute_std_obj(obj_list, history_size);
        }
        if (std_obj <= std_threshold) {
            printf("iter: %d, std_threshold: %.6f\n", iter, std_obj);
            break;
        }

        cur_idx = (x_sol.array() >= 0.5).matrix().cast<float_t>();
        prev_idx = cur_idx;
        cur_obj = compute_cost(prev_idx, _A, _b);

        if (best_bin_obj >= cur_obj) {
            best_bin_obj = cur_obj;
            best_sol = x_sol;
        }
//        norm of x_sol:  60.65611822893488
//        norm of y1:  50.269276
//        norm of y2:  50.269274910227224
//        norm of z1:  336.32866739518795
//        norm of z2:  336.32866739518795
//        norm of cur_idx:  37.52332607858744
        if (printMsg) {
            std::cout << "iter: " << iter << std::endl;
            std::cout << "obj_val: " << obj_val << std::endl;
            std::cout << "cur_obj: " << cur_obj << std::endl;
            std::cout << "norm of x_sol: " << x_sol.norm() << std::endl;
            std::cout << "norm of y1: " << y1.norm() << std::endl;
            std::cout << "norm of y2: " << y2.norm() << std::endl;
            std::cout << "norm of z1: " << z1.norm() << std::endl;
            std::cout << "norm of z2: " << z2.norm() << std::endl;
            std::cout << "norm of cur_idx: " << cur_idx.norm() << std::endl;
            std::cout << "temp1: " << temp1 << std::endl;
            std::cout << "temp2: " << temp2 << std::endl;
            std::cout << "-------------------------------------------------" << std::endl;
        }

#if PRINT_VEC
        std::cout << "x_sol: " << x_sol.format(CommaInitFmt) << std::endl;
            std::cout << "y1: " << y1.format(CommaInitFmt) << std::endl;
            std::cout << "y2: " << y2.format(CommaInitFmt) << std::endl;
            std::cout << "z1: " << z1.format(CommaInitFmt) << std::endl;
            std::cout << "z2: " << z2.format(CommaInitFmt) << std::endl;
#endif
    }
    sol.x_sol = new float_t[n];
    sol.y1 = new float_t[n];
    sol.y2 = new float_t[n];
    sol.best_sol = new float_t[n];

    memcpy(sol.x_sol, x_sol.data(), sizeof(float_t) * n);
    memcpy(sol.y1, y1.data(), sizeof(float_t) * n);
    memcpy(sol.y2, y2.data(), sizeof(float_t) * n);
    memcpy(sol.best_sol, best_sol.data(), sizeof(float_t) * n);
    sol.time_elapsed = time_elapsed;
    return 1;
}

int solver::ADMM_bqp_unconstrained(int matrix_layout, int n, float_t *A, float_t *b, Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto _A       = SparseMatrix(n, n);
    auto _b       = DenseVector(n);

    /* Generating the sparse matrix for A */
    std::vector<Triplet> triplet_list;
    for (int i = 0; i < n * n; i++) {
        if (A[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, A[i]));
    }
    _A.setFromTriplets(triplet_list.begin(), triplet_list.end());
    memcpy(_b.data(), b, n * sizeof(float_t));

    int ret = ADMM_bqp_unconstrained(matrix_layout, n, _A, _b, sol, printMsg);
    return ret;
}

int solver::ADMM_bqp_linear_eq(int matrix_layout, int m, int n, const SparseMatrix &_A, const DenseVector &_b,
                       const SparseMatrix &_C, const DenseVector &_d, Solution& sol, bool printMsg) {
    auto x_sol    = DenseVector(n);
    auto y1       = DenseVector(n);
    auto y2       = DenseVector(n);
    auto z1       = DenseVector(n);
    auto z2       = DenseVector(n);
    auto z3       = DenseVector(m);
    auto prev_idx = DenseVector(n);
    auto best_sol = DenseVector(n);
    auto temp_vec = DenseVector(n);
    auto _2A      = SparseMatrix(n, n);
    auto temp_mat = SparseMatrix(n, n);
    auto cur_idx  = DenseVector(n);
    auto Csq      = SparseMatrix(n, n);
    float_t cur_obj;
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteLUT<float_t>> cg;
    cg.setMaxIterations(pcg_maxiters);
    cg.setTolerance(pcg_tol);

    _2A = 2 * _A;

    Csq = _C.transpose() * _C;

    if (x0 != nullptr) {
        memcpy(x_sol.data(), x0, sizeof(float_t) * n);
    } else {
        memset(x_sol.data(), 0, sizeof(float_t) * n);
    }
    y1 = x_sol;
    y2 = x_sol;
    prev_idx = x_sol;
    best_sol = x_sol;

    float_t rho1 = initial_rho;
    float_t rho2 = initial_rho;
    float_t rho3 = initial_rho;
    std::vector<float_t> obj_list; /* Stores the objective value calculated during each iteration */
    float_t std_obj = 1;

    float_t best_bin_obj = compute_cost(x_sol, _A, _b);

    long time_elapsed = 0;
    for (int iter = 0; iter < max_iters; iter++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        temp_vec = x_sol + z1 / rho1;

        /* Project vector on [0, 1] box */
        project_box(n, temp_vec, y1);

        temp_vec = x_sol + z2 / rho2;

        /* Project vector on shifted lp box */
        project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

        temp_mat = _2A + rho3 * Csq;
        for (int i = 0; i < n; i++) {
            temp_mat.coeffRef(i, i) += rho1 + rho2;
        }
        temp_vec = -(_b + z1 + z2 + _C.transpose() * z3) + rho1 * y1 + rho2 * y2 + rho3 * _C.transpose() * _d;

        /* The solveWithGuess function does not seems to the the tolerance condition in the beginning */
        if ((temp_mat * y1 - temp_vec).norm() < cg.tolerance() * _b.norm()) {
            x_sol = y1;
        } else {
            cg.compute(temp_mat);
            x_sol = cg.solveWithGuess(temp_vec, y1);
        }


        z1 = z1 + gamma_val * rho1 * (x_sol - y1);
        z2 = z2 + gamma_val * rho2 * (x_sol - y2);
        z3 = z3 + gamma_val * rho3 * (_C * x_sol - _d);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time_elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

        if ((iter+1) % rho_change_step == 0) {
            rho1 = learning_fact * rho1;
            rho2 = learning_fact * rho2;
            rho3 = learning_fact * rho3;
            gamma_val = std::max(gamma_val * gamma_factor, float_t(1));
        }

        float_t temp0 = std::max(x_sol.norm(), float_t(2.2204e-16));
        float_t temp1 = (x_sol-y1).norm() / temp0;
        float_t temp2 = (x_sol-y2).norm() / temp0;
        if (temp1 <= stop_threshold && temp2 <= stop_threshold) {
            printf("iter: %d, stop_threshold: %.6f", iter, std::max(temp1,temp2));
            break;
        }
        float_t obj_val = compute_cost(x_sol,_A,_b);
        obj_list.push_back(obj_val);
        if (obj_list.size() >= history_size) {
            std_obj = compute_std_obj(obj_list, history_size);
        }

        if (std_obj <= std_threshold) {
            printf("iter: %d, std_threshold: %.6f", iter, std_obj);
            break;
        }

        cur_idx = (x_sol.array() >= 0.5).matrix().cast<float_t>();
        prev_idx = cur_idx;
        cur_obj = compute_cost(prev_idx, _A, _b);

        if (best_bin_obj >= cur_obj) {
            best_bin_obj = cur_obj;
            best_sol = x_sol;
        }

        if (printMsg) {
            std::cout << "iter: " << iter << std::endl;
            std::cout << "obj_val: " << obj_val << std::endl;
            std::cout << "cur_obj: " << cur_obj << std::endl;
            std::cout << "temp1: " << temp1 << std::endl;
            std::cout << "temp2: " << temp2 << std::endl;
            std::cout << "-------------------------------------------------" << std::endl;
        }

#if PRINT_VEC
        std::cout << "x_sol: " << x_sol.format(CommaInitFmt) << std::endl;
            std::cout << "y1: " << y1.format(CommaInitFmt) << std::endl;
            std::cout << "y2: " << y2.format(CommaInitFmt) << std::endl;
            std::cout << "z1: " << z1.format(CommaInitFmt) << std::endl;
            std::cout << "z2: " << z2.format(CommaInitFmt) << std::endl;
#endif
    }

    sol.x_sol = new float_t[n];
    sol.y1 = new float_t[n];
    sol.y2 = new float_t[n];
    sol.best_sol = new float_t[n];

    memcpy(sol.x_sol, x_sol.data(), sizeof(float_t) * n);
    memcpy(sol.y1, y1.data(), sizeof(float_t) * n);
    memcpy(sol.y2, y2.data(), sizeof(float_t) * n);
    memcpy(sol.best_sol, best_sol.data(), sizeof(float_t) * n);
    sol.time_elapsed = time_elapsed;
    return 1;

}

int solver::ADMM_bqp_linear_eq(int matrix_layout, int n, int m, float_t *A, float_t *b, float_t *C, float_t *d,
                               Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto _A       = SparseMatrix(n, n);
    auto _b       = DenseVector(n);
    auto _C       = SparseMatrix(m, n);
    auto _d       = DenseVector(m);

    /* Generating the sparse matrix for A */
    std::vector<Triplet> triplet_list; /* Maybe one can reserve some space here */
    for (int i = 0; i < n * n; i++) {
        if (A[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, A[i]));
    }
    _A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    /* Generating the sparse matrix for C */
    triplet_list.clear();
    for (int i = 0; i < m * n; i++) {
        if (C[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, C[i]));
    }
    _C.setFromTriplets(triplet_list.begin(), triplet_list.end());

    memcpy(_b.data(), b, n * sizeof(float_t));
    memcpy(_d.data(), d, m * sizeof(float_t));

    int ret = ADMM_bqp_linear_eq(matrix_layout, n, m, _A, _b, _C, _d, sol, printMsg);
    return ret;
}

int solver::ADMM_bqp_linear_ineq(int matrix_layout, int n, int m, const SparseMatrix &_A, const DenseVector &_b,
                         const SparseMatrix &_E, const DenseVector &_f, Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto x_sol    = DenseVector(n);
    auto y1       = DenseVector(n);
    auto y2       = DenseVector(n);
    auto y3       = DenseVector(m);
    auto z1       = DenseVector(n);
    auto z2       = DenseVector(n);
    auto z4       = DenseVector(m);
    auto prev_idx = DenseVector(n);
    auto best_sol = DenseVector(n);
    auto temp_vec = DenseVector(n);
    auto _2A      = SparseMatrix(n, n);
    auto temp_mat = SparseMatrix(n, n);
    auto cur_idx  = DenseVector(n);
    auto Esq      = SparseMatrix(n, n);
    float_t cur_obj;
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteLUT<float_t>> cg;
    cg.setMaxIterations(pcg_maxiters);
    cg.setTolerance(pcg_tol);

    _2A = 2 * _A;

    Esq = _E.transpose() * _E;

    if (x0 != nullptr) {
        memcpy(x_sol.data(), x0, sizeof(float_t) * n);
    } else {
        memset(x_sol.data(), 0, sizeof(float_t) * n);
    }
    y1 = x_sol;
    y2 = x_sol;
    y3 = _f - _E * x_sol;
    prev_idx = x_sol;
    best_sol = x_sol;

    float_t rho1 = initial_rho;
    float_t rho2 = initial_rho;
    float_t rho4 = initial_rho;
    std::vector<float_t> obj_list; /* Stores the objective value calculated during each iteration */
    float_t std_obj = 1;

    float_t best_bin_obj = compute_cost(x_sol, _A, _b);
    std::cout << "best_bin_obj: " << best_bin_obj << std::endl;
    cur_obj = best_bin_obj;

    long time_elapsed = 0;
    for (int iter = 0; iter < max_iters; iter++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        temp_vec = x_sol + z1 / rho1;

        /* Project vector on [0, 1] box */
        project_box(n, temp_vec, y1);

        temp_vec = x_sol + z2 / rho2;

        /* Project vector on shifted lp box */
        project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

        y3 = (_f - _E * x_sol - z4 / rho4).array().max(0).matrix();

        temp_mat = _2A + rho4 * Esq;
        for (int i = 0; i < n; i++) {
            temp_mat.coeffRef(i, i) += rho1 + rho2;
        }
        temp_vec = -(_b + z1 + z2 + _E.transpose() * z4) + rho1 * y1 + rho2 * y2 + rho4 * _E.transpose() * (_f - y3);

        /* The solveWithGuess function does not seems to the the tolerance condition in the beginning */
        if ((temp_mat * y1 - temp_vec).norm() < cg.tolerance() * _b.norm()) {
            x_sol = y1;
        } else {
            cg.compute(temp_mat);
            x_sol = cg.solveWithGuess(temp_vec, y1);
        }


        z1 = z1 + gamma_val * rho1 * (x_sol - y1);
        z2 = z2 + gamma_val * rho2 * (x_sol - y2);
        z4 = z4 + gamma_val * rho4 * (_E * x_sol + y3 - _f);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time_elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

        if ((iter+2) % rho_change_step == 0) {
            rho1 = learning_fact * rho1;
            rho2 = learning_fact * rho2;
            rho4 = learning_fact * rho4;
            gamma_val = std::max(gamma_val * gamma_factor, float_t(1));
        }

        float_t temp0 = std::max(x_sol.norm(), float_t(2.2204e-16));
        float_t temp1 = (x_sol-y1).norm() / temp0;
        float_t temp2 = (x_sol-y2).norm() / temp0;
        if (temp1 <= stop_threshold && temp2 <= stop_threshold) {
            printf("iter: %d, stop_threshold: %.6f", iter, std::max(temp1,temp2));
            break;
        }
        float_t obj_val = compute_cost(x_sol,_A,_b);
        obj_list.push_back(obj_val);
        if (obj_list.size() >= history_size) {
            std_obj = compute_std_obj(obj_list, history_size);
        }
        if (std_obj <= std_threshold) {
            printf("iter: %d, std_threshold: %.6f", iter, std_obj);
            break;
        }

        cur_idx = (x_sol.array() >= 0.5).matrix().cast<float_t>();
        prev_idx = cur_idx;
        cur_obj = compute_cost(prev_idx, _A, _b);

        if (best_bin_obj >= cur_obj) {
            best_bin_obj = cur_obj;
            best_sol = x_sol;
        }

        std::cout << "iter: " << iter << std::endl;
        if (printMsg) {
            std::cout << "obj_val: " << obj_val << std::endl;
            std::cout << "cur_obj: " << cur_obj << std::endl;
            std::cout << "temp1: " << temp1 << std::endl;
            std::cout << "temp2: " << temp2 << std::endl;
        }
#if PRINT_VEC
        std::cout << "x_sol: " << x_sol.format(CommaInitFmt) << std::endl;
            std::cout << "y1: " << y1.format(CommaInitFmt) << std::endl;
            std::cout << "y2: " << y2.format(CommaInitFmt) << std::endl;
            std::cout << "z1: " << z1.format(CommaInitFmt) << std::endl;
            std::cout << "z2: " << z2.format(CommaInitFmt) << std::endl;
#endif
        std::cout << "-------------------------------------------------" << std::endl;
    }
    sol.x_sol = new float_t[n];
    sol.y1 = new float_t[n];
    sol.y2 = new float_t[n];
    sol.best_sol = new float_t[n];

    memcpy(sol.x_sol, x_sol.data(), sizeof(float_t) * n);
    memcpy(sol.y1, y1.data(), sizeof(float_t) * n);
    memcpy(sol.y2, y2.data(), sizeof(float_t) * n);
    memcpy(sol.best_sol, best_sol.data(), sizeof(float_t) * n);
    sol.time_elapsed = time_elapsed;
    return 1;
}

int solver::ADMM_bqp_linear_ineq(int matrix_layout, int n, int m, float_t *A, float_t *b, float_t *E, float_t *f,
                                 Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto _A       = SparseMatrix(n, n);
    auto _b       = DenseVector(n);
    auto _E       = SparseMatrix(m, n);
    auto _f       = DenseVector(m);

    /* Generating the sparse matrix for A */
    std::vector<Triplet> triplet_list; /* Maybe one can reserve some space here */
    for (int i = 0; i < n * n; i++) {
        if (A[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, A[i]));
    }
    _A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    /* Generating the sparse matrix for C */
    triplet_list.clear();
    for (int i = 0; i < m * n; i++) {
        if (E[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, E[i]));
    }
    _E.setFromTriplets(triplet_list.begin(), triplet_list.end());

    memcpy(_b.data(), b, n * sizeof(float_t));
    memcpy(_f.data(), f, m * sizeof(float_t));

    int ret = ADMM_bqp_linear_ineq(matrix_layout, n, m, _A, _b, _E, _f, sol, printMsg);
    return ret;
}

int solver::ADMM_bqp_linear_eq_and_uneq(int matrix_layout, int l, int m, int n, const SparseMatrix &_A, const DenseVector &_b,
                                const SparseMatrix &_C, const DenseVector &_d, const SparseMatrix &_E,
                                const DenseVector &_f, Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto x_sol    = DenseVector(n);
    auto y1       = DenseVector(n);
    auto y2       = DenseVector(n);
    auto y3       = DenseVector(m);
    auto z1       = DenseVector(n);
    auto z2       = DenseVector(n);
    auto z3       = DenseVector(l);
    auto z4       = DenseVector(m);
    auto prev_idx = DenseVector(n);
    auto best_sol = DenseVector(n);
    auto temp_vec = DenseVector(n);
    auto _2A      = SparseMatrix(n, n);
    auto temp_mat = SparseMatrix(n, n);
    auto Csq      = SparseMatrix(n, n);
    auto cur_idx  = DenseVector(n);
    auto Esq      = SparseMatrix(n, n);
    float_t cur_obj;
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteLUT<float_t>> cg;
    cg.setMaxIterations(pcg_maxiters);
    cg.setTolerance(pcg_tol);

    _2A = 2 * _A;
    Esq = _E.transpose() * _E;
    Csq = _C.transpose() * _C;

    if (x0 != nullptr) {
        memcpy(x_sol.data(), x0, sizeof(float_t) * n);
    } else {
        memset(x_sol.data(), 0, sizeof(float_t) * n);
    }
    y1 = x_sol;
    y2 = x_sol;
    y3 = _f - _E * x_sol;
    prev_idx = x_sol;
    best_sol = x_sol;

    float_t rho1 = initial_rho;
    float_t rho2 = initial_rho;
    float_t rho3 = initial_rho;
    float_t rho4 = initial_rho;
    std::vector<float_t> obj_list; /* Stores the objective value calculated during each iteration */
    float_t std_obj = 1;

    float_t best_bin_obj = compute_cost(x_sol, _A, _b);

    long time_elapsed = 0;
    for (int iter = 0; iter < max_iters; iter++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        temp_vec = x_sol + z1 / rho1;

        /* Project vector on [0, 1] box */
        project_box(n, temp_vec, y1);

        temp_vec = x_sol + z2 / rho2;

        /* Project vector on shifted lp box */
        project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

        y3 = (_f - _E * x_sol - z4 / rho4).array().max(0).matrix();

        temp_mat = _2A + rho3 * Csq + rho4 * Esq;
        for (int i = 0; i < n; i++) {
            temp_mat.coeffRef(i, i) += rho1 + rho2;
        }
        temp_vec = -(_b + z1 + z2 + _C.transpose() * z3 + _E.transpose() * z4) + rho1 * y1 + rho2 * y2 + rho3 * _C.transpose() * _d
                   + rho4 * _E.transpose() * (_f - y3);

        /* The solveWithGuess function does not seems to the the tolerance condition in the beginning */
        if ((temp_mat * y1 - temp_vec).norm() < cg.tolerance() * _b.norm()) {
            x_sol = y1;
        } else {
            cg.compute(temp_mat);
            x_sol = cg.solveWithGuess(temp_vec, y1);
        }


        z1 = z1 + gamma_val * rho1 * (x_sol - y1);
        z2 = z2 + gamma_val * rho2 * (x_sol - y2);
        z3 = z3 + gamma_val * rho3 * (_C * x_sol - _d);
        z4 = z4 + gamma_val * rho4 * (_E * x_sol + y3 - _f);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        time_elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

        if ((iter+2) % rho_change_step == 0) {
            rho1 = learning_fact * rho1;
            rho2 = learning_fact * rho2;
            rho4 = learning_fact * rho4;
            gamma_val = std::max(gamma_val * gamma_factor, float_t(1));
        }

        float_t temp0 = std::max(x_sol.norm(), float_t(2.2204e-16));
        float_t temp1 = (x_sol-y1).norm() / temp0;
        float_t temp2 = (x_sol-y2).norm() / temp0;
        if (temp1 <= stop_threshold && temp2 <= stop_threshold) {
            printf("iter: %d, stop_threshold: %.6f", iter, std::max(temp1,temp2));
            break;
        }

        float_t obj_val = compute_cost(x_sol,_A,_b);
        obj_list.push_back(obj_val);
        if (obj_list.size() >= history_size) {
            std_obj = compute_std_obj(obj_list, history_size);
        }
        if (std_obj <= std_threshold) {
            printf("iter: %d, std_threshold: %.6f", iter, std_obj);
            break;
        }

        cur_idx = (x_sol.array() >= 0.5).matrix().cast<float_t>();
        prev_idx = cur_idx;
        cur_obj = compute_cost(prev_idx, _A, _b);

        if (best_bin_obj >= cur_obj) {
            best_bin_obj = cur_obj;
            best_sol = x_sol;
        }
        std::cout << "iter: " << iter << std::endl;
        if (printMsg) {
            std::cout << "obj_val: " << obj_val << std::endl;
            std::cout << "cur_obj: " << cur_obj << std::endl;
            std::cout << "temp1: " << temp1 << std::endl;
            std::cout << "temp2: " << temp2 << std::endl;
        }
#if PRINT_VEC
        std::cout << "x_sol: " << x_sol.format(CommaInitFmt) << std::endl;
            std::cout << "y1: " << y1.format(CommaInitFmt) << std::endl;
            std::cout << "y2: " << y2.format(CommaInitFmt) << std::endl;
            std::cout << "z1: " << z1.format(CommaInitFmt) << std::endl;
            std::cout << "z2: " << z2.format(CommaInitFmt) << std::endl;
#endif
        std::cout << "-------------------------------------------------" << std::endl;
    }
    sol.x_sol = new float_t[n];
    sol.y1 = new float_t[n];
    sol.y2 = new float_t[n];
    sol.best_sol = new float_t[n];

    memcpy(sol.x_sol, x_sol.data(), sizeof(float_t) * n);
    memcpy(sol.y1, y1.data(), sizeof(float_t) * n);
    memcpy(sol.y2, y2.data(), sizeof(float_t) * n);
    memcpy(sol.best_sol, best_sol.data(), sizeof(float_t) * n);
    sol.time_elapsed = time_elapsed;
    return 1;

}

int solver::ADMM_bqp_linear_eq_and_uneq(int matrix_layout, int l, int m, int n, float_t *A, float_t *b, float_t *C, float_t *d,
                                        float_t *E, float_t *f, Solution& sol, bool printMsg) {
    /* Initialize the parameters */
    auto _A       = SparseMatrix(n, n);
    auto _b       = DenseVector(n);
    auto _C       = SparseMatrix(l, n);
    auto _d       = DenseVector(l);
    auto _E       = SparseMatrix(m, n);
    auto _f       = DenseVector(m);

    /* Generating the sparse matrix for A */
    std::vector<Triplet> triplet_list; /* Maybe one can reserve some space here */
    for (int i = 0; i < n * n; i++) {
        if (A[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, A[i]));
    }
    _A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    /* Generating the sparse matrix for C */
    triplet_list.clear();
    for (int i = 0; i < m * n; i++) {
        if (E[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, E[i]));
    }
    _E.setFromTriplets(triplet_list.begin(), triplet_list.end());

    /* Generating the sparse matrix for C */
    triplet_list.clear();
    for (int i = 0; i < m * n; i++) {
        if (C[i] == 0.0) {
            continue;
        }
        int row = i % n;
        int col = i / n;
        triplet_list.push_back(Triplet(row, col, C[i]));
    }
    _C.setFromTriplets(triplet_list.begin(), triplet_list.end());

    memcpy(_b.data(), b, n * sizeof(float_t));
    memcpy(_d.data(), d, l * sizeof(float_t));
    memcpy(_f.data(), f, m * sizeof(float_t));

    int ret = ADMM_bqp_linear_eq_and_uneq(matrix_layout, l, m, n, _A, _b, _C, _d, _E, _f, sol, printMsg);
    return ret;
}
