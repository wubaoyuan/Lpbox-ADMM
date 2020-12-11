//
// Created by squintingsmile on 2020/9/22.
//

#ifndef CPP_SOLVER_H
#define CPP_SOLVER_H

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define float_t double
#define Dynamic Eigen::Dynamic
#define PRINT_VEC false
typedef Eigen::SparseMatrix<float_t> SparseMatrix;             /* Sparse matrix type */
typedef Eigen::Matrix<float_t, Dynamic, Dynamic> DenseMatrix;  /* Dense matrix type */
typedef Eigen::Triplet<float_t> Triplet;                       /* Triplet type */
typedef Eigen::Matrix<float_t, Dynamic, 1> DenseVector;        /* Dense vector type */

enum Preconditioner {
    INCOMPLETE_CHOLESKY = 0
};

struct Solution {
    float_t *best_sol;
    float_t *x_sol;
    float_t *y1;
    float_t *y2;
    double time_elapsed;
};

/* Returns c = a + b */
inline void vec_add(int n, double *a, double *b, double *c) {
    for (int i = 0; i < n; i++) {
        c[i] = a[i] + b[i];
    }
}

/* Returns c = a - b */
inline void vec_sub(int n, double *a, double *b, double *c) {
    for (int i = 0; i < n; i++) {
        c[i] = a[i] - b[i];
    }
}

/* Set all the elements of a vector a to certain value */
inline void vec_set(int n, double* a, const double val) {
    for (int i = 0; i < n; i++) {
        a[i] = val;
    }
}

class solver {
private:
    float_t std_dev(std::vector<float_t> arr, size_t begin, size_t end); /* Calculate the standard deviation */

    float_t stop_threshold;
    float_t std_threshold;
    int     max_iters;
    float_t initial_rho;
    int     rho_change_step;
    float_t gamma_val;
    float_t learning_fact;
    float_t history_size;
    float_t projection_lp;
    float_t gamma_factor;
    float_t pcg_tol;
    float_t rho_upper_limit;
    float_t rel_tol;
    int     pcg_maxiters;
    float_t *x0 = nullptr;

    void project_box(int n, const double *x, double *y); /* Project the vector x on the box [0,1]^n and
                                                            * return the projection in y */

    void project_box(int n, const DenseVector &x, DenseVector &y);

    void project_shifted_Lp_ball(int n, double *x, int p, double *y); /* Project the vector x on the shifted
                                                                         * norm box and return the projection in y */

    void project_shifted_Lp_ball(int n, const DenseVector &x, int p, DenseVector &y);

    double compute_cost(int m, int n, double *x, double *A, double *b); /* x^TAx + b^Tx */

    float_t compute_cost(const DenseVector &x, const SparseMatrix &A, const DenseVector &b);

    float_t compute_std_obj(std::vector<float_t> obj_list, int history_size); /* Compute the standard deviation of the objective
                                                                               * values of the several latest iterations */

    int conjugate_gradient_with_guess(const SparseMatrix& A, const DenseVector& b, const DenseVector& init_guess,
                                      DenseVector &result, double tol, int max_iter, bool printMsg=false);

    int conjugate_gradient_with_guess_pre(const SparseMatrix &A, const DenseVector &b, const DenseVector &init_guess,
                                              DenseVector &result, double tol, int max_iter,
                                              const Preconditioner preconditioner, bool printMsg=false);

public:

    void ADMM_bqp_unconstrained_init();

    void ADMM_bqp_linear_eq_init();

    void ADMM_bqp_linear_ineq_init();

    void ADMM_bqp_linear_eq_and_uneq_init();

    int ADMM_bqp_unconstrained(int matrix_layout, int n, const SparseMatrix &_A, const DenseVector &_b, Solution& sol, bool printMsg);

    int ADMM_bqp_unconstrained(int matrix_layout, int n, float_t *A, float_t *b, Solution& sol, bool printMsg);

    int ADMM_bqp_linear_eq(int matrix_layout, int m, int n, const SparseMatrix &_A, const DenseVector &_b,
                           const SparseMatrix &_C, const DenseVector &_d, Solution& sol, bool printMsg);

    int ADMM_bqp_linear_eq(int matrix_layout, int m, int n, float_t *A, float_t *b, float_t *C, float_t *d, Solution& sol, bool printMsg);

    int ADMM_bqp_linear_ineq(int matrix_layout, int n, int m, const SparseMatrix &_A, const DenseVector &_b,
                             const SparseMatrix &_E, const DenseVector &_f, Solution& sol, bool printMsg);

    int ADMM_bqp_linear_ineq(int matrix_layout, int n, int m, float_t *A, float_t *b, float_t *E, float_t *f,
                             Solution& sol, bool printMsg);

    int ADMM_bqp_linear_eq_and_uneq(int matrix_layout, int l, int m, int n, const SparseMatrix &_A, const DenseVector &_b,
                                    const SparseMatrix &_C, const DenseVector &_d, const SparseMatrix &_E,
                                    const DenseVector &_f, Solution& sol, bool printMsg);

    int ADMM_bqp_linear_eq_and_uneq(int matrix_layout, int l, int m, int n, float_t *A, float_t *b, float_t *C, float_t *d,
                                    float_t *E, float_t *f, Solution& sol, bool printMsg);

    void set_stop_threshold(float_t stop_threshold) {
        this->stop_threshold = stop_threshold;
    }

    void set_std_threshold(float_t std_threshold) {
        this->std_threshold = std_threshold;
    }

    void set_max_iters(int max_iters) {
        this->max_iters = max_iters;
    }

    void set_initial_rho(float_t initial_rho) {
        this->initial_rho = initial_rho;
    }

    void set_rho_change_step(float_t rho_change_step) {
        this->rho_change_step = rho_change_step;
    }

    void set_gamma_val(float_t gamma_val) {
        this->gamma_val = gamma_val;
    }

    void set_learning_fact(float_t learning_fact) {
        this->learning_fact = learning_fact;
    }

    void set_history_size(float_t history_size) {
        this->history_size = history_size;
    }

    void set_projection_lp(float_t projection_lp) {
        this->projection_lp = projection_lp;
    }

    void set_gamma_factor(float_t gamma_factor) {
        this->gamma_factor = gamma_factor;
    }

    void set_pcg_tol(float_t pcg_tol) {
        this->pcg_tol = pcg_tol;
    }

    void set_rho_upper_limit(float_t rho_upper_limit) {
        this->rho_upper_limit = rho_upper_limit;
    }

    void set_rel_tol(float_t rel_tol) {
        this->rel_tol = rel_tol;
    }

    void set_pcg_maxiters(int pcg_maxiters) {
        this->pcg_maxiters = pcg_maxiters;
    }

    void set_x0(float_t *x0) {
        this->x0 = x0;
    }

};

#endif //CPP_SOLVER_H
