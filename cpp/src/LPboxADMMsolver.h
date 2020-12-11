//
// Created by squintingsmile on 2020/9/22.
//

#ifndef CPP_LPBOXADMMSOLVER_H
#define CPP_LPBOXADMMSOLVER_H

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define Dynamic Eigen::Dynamic
#define PRINT_VEC false
#define PRINT_SHAPE(mat) printf("%s col: %ld, row: %ld\n", #mat, mat.cols(), mat.rows())
#define PRINT_NORM(mat) printf("%s col: %ld, row: %ld, norm: %lf\n", #mat, mat.cols(), mat.rows(), mat.norm())
typedef double _double_t;
typedef Eigen::SparseMatrix<_double_t, Eigen::RowMajorBit> SparseMatrix;             /* Sparse matrix type. The default
																					  *	is row major to speed up matrix
																					  *	multiplication */
typedef Eigen::Matrix<_double_t, Dynamic, Dynamic> DenseMatrix;  /* Dense matrix type */
typedef Eigen::Triplet<_double_t> Triplet;                       /* Triplet type */
typedef Eigen::Matrix<_double_t, Dynamic, 1> DenseVector;        /* Dense vector type */
typedef Eigen::Matrix<int, Dynamic, Dynamic> DenseIntMatrix;     /* Matrix with integer values */
typedef Eigen::Matrix<int, Dynamic, 1> DenseIntVector;			 /* Vector with integer values */

enum PreconditionerType {
	DIAGONAL = 0,
	INCOMPLETE_CHOLESKY = 1,
	IDENTITY = 2
};

/* Stores the solution return by the algorithm */
struct Solution {
	DenseVector *best_sol;
	DenseVector *x_sol;
	DenseVector *y1;
	DenseVector *y2;
	long time_elapsed;
};

enum problem_t {
    unconstrained,
    equality,
    inequality,
    equality_and_inequality
};

/* Storing the pointers to the matrices for matrix_expression */
struct MatrixInfo {
    problem_t problem_type;
    int n;
    const SparseMatrix *A = nullptr;
    const DenseVector *b = nullptr;
    const DenseVector *x0 = nullptr;
    int m;
    const SparseMatrix *C = nullptr;
    const DenseVector *d = nullptr;
    int l;
    const SparseMatrix *E = nullptr;
    const DenseVector *f = nullptr;
};

/* Storing the instruction for the updating behavior of the algorithm */
struct SolverInstruction{
    problem_t problem_type;
    int update_y3;
    int update_z3;
    int update_z4;
	int update_rho3;
	int update_rho4;
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

class LPboxADMMsolver {
	private:

		/* The stopping threashold for the convergence condition 
		 * (see .cpp file for termination condition)
		 */
		_double_t stop_threshold;

		/* If the standard deviation of objective values within the 
		 * history size is less than the std_threashold, the algorithm
		 * terminates
		 */
		_double_t std_threshold;

		/* Maximum iteration number of the algortithm */
		int max_iters;
		_double_t initial_rho;

		/* This defines how many iterations will it take before the rho value is updated */
		int rho_change_step;
		_double_t gamma_val;
		_double_t learning_fact;

		/* The number of previous objective value that are used to calculate standard deviation */
		_double_t history_size;

		/* The p-norm that the algorithm uses to update y2. Default is 2-norm */
		_double_t projection_lp;
		_double_t gamma_factor;

		/* The tolerance of conjugate gradient */
		_double_t pcg_tol;
		int pcg_maxiters;

		/* The upper limit for rho increasing */
		_double_t rho_upper_limit;
		_double_t rel_tol;

		/* Determine whether the solver log the intermediate steps */
		int does_log = 0;
		std::string log_file_path;

		/* The preconditioner used by the conjugate gradient methods (currently unused) */
		PreconditionerType preconditioner;

		/* Calculate the standard deviation */
		_double_t std_dev(std::vector<_double_t>& arr, size_t begin, size_t end); 

		/* Project the vector value greater than `greater_val` to `set_val` */
		void project_vec_greater_than(DenseVector &x, DenseVector &res, const int &greater_val, const int &set_val);

		/* Project the vector value less than `less_val` to `set_val` */
		void project_vec_less_than(DenseVector &x, DenseVector &res, const int &less_val, const int &set_val);

		/* Project the vector x on the box [0,1]^n and return the projection in y */
		void project_box(int n, const double *x, double *y);

		void project_box(int n, const DenseVector &x, DenseVector &y);

		/* Project the vector x on the shifted norm box and return the projection in y */
		void project_shifted_Lp_ball(int n, double *x, int p, double *y);

		/* Project the vector x on the shifted norm box [0, 1]^n and return the projection in y */
		void project_shifted_Lp_ball(int n, const DenseVector &x, int p, DenseVector &y);

		/* Compute the cost x^TAx + b^Tx*/
		double compute_cost(int m, int n, double *x, double *A, double *b);

		_double_t compute_cost(const DenseVector &x, const SparseMatrix &A, const DenseVector &b);
		
		_double_t compute_cost(const DenseVector &x, const SparseMatrix &A, const DenseVector &b, 
				DenseVector &temp_vec_for_mat_mul);

		/* Compute the standard deviation of the objective values of the several latest iterations */
		_double_t compute_std_obj(std::vector<_double_t> obj_list, int history_size);	

	public:

		/* Default parameter initialization for unconstraint problem */
		void ADMM_bqp_unconstrained_init();

		void ADMM_bqp_linear_eq_init();

		void ADMM_bqp_linear_ineq_init();

		void ADMM_bqp_linear_eq_and_uneq_init();

		/**
		 * This function solves the following optimization problem 
		 *			min_x (x^T A x + b^T x) 
		 * such that x is {0,1}^n ADMM update steps with x0 is feasible and binary.
		 * @param matrix_layout: The layout of the matrix (currently not used)
		 * @param n:  The length of the b vector and the edges of matrix A
		 * @param _A: The matrix A
		 * @param _b: The vector b
		 * @param sol: Stores the returned solution values
		 * @param printMsg: If this is set to true, print the intermidiate messages
		 */	
		int ADMM_bqp_unconstrained(int n, const SparseMatrix &_A, const DenseVector &_b, 
				const DenseVector &x0, Solution& sol);

		/* The legacy version originally developed. The experiment shows that it is actually faster */
		int ADMM_bqp_unconstrained_legacy(int n, const SparseMatrix &_A, const DenseVector &_b, 
				const DenseVector &x0, Solution& sol);

		/* This provides an interface to use double array as the algorithm input, the matrix A here is a column major
		 * dense matrix */
		int ADMM_bqp_unconstrained(int n, _double_t *A, _double_t *b, _double_t *x0, Solution& sol);

		/**
		 * This function solves the following optimization problem 
		 *			min_x (x^T A x + b^T x)
		 *				s.t.	Cx = d
		 * such that x is {0,1}^n ADMM update steps with x0 is feasible and binary.
		 * @param matrix_layout: The layout of the matrix (currently not used)
		 * @param n:  The length of the b vector and the edges of matrix A
		 * @param _A: The matrix A
		 * @param _b: The vector b
		 * @param m:  The length of the d vector and the edges of matrix C
		 * @param _C: The matrix C
		 * @param _d: The matrix d
		 * @param sol: Stores the returned solution values
		 * @param printMsg: If this is set to true, print the intermidiate messages
		 */	
		int ADMM_bqp_linear_eq(int n, const SparseMatrix &_A, const DenseVector &_b, const DenseVector &x0,
				int m, const SparseMatrix &_C, const DenseVector &_d, Solution& sol);

		int ADMM_bqp_linear_eq(int n, _double_t *A, _double_t *b, _double_t *x0, int m, _double_t *C, 
				_double_t *d, Solution& sol);

		/* This solves for the case when the matrix A is the kronckner product of kron and _A. This assumes that the 
		 * kron matrix does not have entries that are on the same row or the same column. For example, kron can be
		 * [ 0, 1;  [ 1, 0;	    but can't be [ 0, 1;
		 *	 1, 0 ]   0, 1 ]				   0, 1 ]
		 * If kron is [ 1, 0;  , then the resulting matrix is [ _A, 0;
		 *				0, 1 ]									0, _A ]
		 * 
		 * TODO: Implement this function. Note that there is already workaround method to solve the matrix by matrix
		 * multiplication function so this function might not be needed anymore
		 */
		int ADMM_bqp_linear_eq_kron(int n, const SparseMatrix &_A, const DenseMatrix &kron, const DenseVector &_b,
				int m, const SparseMatrix &_C, const DenseVector &_d, Solution& sol);

		/**
		 * This function solves the following optimization problem 
		 *			min_x (x^T A x + b^T x)
		 *				s.t.	Ex <= f
		 * such that x is {0,1}^n ADMM update steps with x0 is feasible and binary.
		 * @param matrix_layout: The layout of the matrix (currently not used)
		 * @param n:  The length of the b vector and the edges of matrix A
		 * @param _A: The matrix A
		 * @param _b: The vector b
		 * @param l:  The length of the f vector and the edges of matrix E
		 * @param _E: The matrix E
		 * @param _f: The vector f
		 * @param sol: Stores the returned solution values
		 * @param printMsg: If this is set to true, print the intermidiate messages
		 */	
		int ADMM_bqp_linear_ineq(int n, const SparseMatrix &_A, const DenseVector &_b, const DenseVector &x0,
				int l, const SparseMatrix &_E, const DenseVector &_f, Solution& sol);

		int ADMM_bqp_linear_ineq(int n, _double_t *A, _double_t *b, _double_t *x0, int l, _double_t *E, _double_t *f,
				Solution& sol);

		/**
		 * This function solves the following optimization problem 
		 *			min_x (x^T A x + b^T x)
		 *				s.t.	Cx = d
		 *						Ex <= f
		 * such that x is {0,1}^n ADMM update steps with x0 is feasible and binary.
		 * @param matrix_layout: The layout of the matrix (currently not used)
		 * @param n:  The length of the b vector and the edges of matrix A
		 * @param _A: The matrix A
		 * @param _b: The vector b
		 * @param m:  The length of the d vector and the edges of matrix C
		 * @param _C: The matrix C
		 * @param _d: The vector d
		 * @param l:  The length of the f vector and the edges of matrix E
		 * @param _E: The matrix E
		 * @param _f: The vector f
		 * @param sol: Stores the returned solution values
		 * @param printMsg: If this is set to true, print the intermidiate messages
		 */	
		int ADMM_bqp_linear_eq_and_uneq(int n, const SparseMatrix &_A, const DenseVector &_b, 
				const DenseVector &x0, int m, const SparseMatrix &_C, const DenseVector &_d, int l, const SparseMatrix &_E,
				const DenseVector &_f, Solution& sol);

		int ADMM_bqp_linear_eq_and_uneq(int n, _double_t *A, _double_t *b, _double_t *x0, 
				int m, _double_t *C, _double_t *d, int l, _double_t *E, _double_t *f, Solution& sol);
		
		/* The actual implementation of the BQP function */
		int ADMM_bqp(MatrixInfo matrix_info, SolverInstruction solver_instruction, Solution& sol);


		void set_stop_threshold(_double_t stop_threshold) {
			this->stop_threshold = stop_threshold;
		}

		void set_std_threshold(_double_t std_threshold) {
			this->std_threshold = std_threshold;
		}

		void set_max_iters(int max_iters) {
			this->max_iters = max_iters;
		}

		void set_initial_rho(_double_t initial_rho) {
			this->initial_rho = initial_rho;
		}

		void set_rho_change_step(_double_t rho_change_step) {
			this->rho_change_step = rho_change_step;
		}

		void set_gamma_val(_double_t gamma_val) {
			this->gamma_val = gamma_val;
		}

		void set_learning_fact(_double_t learning_fact) {
			this->learning_fact = learning_fact;
		}

		void set_history_size(_double_t history_size) {
			this->history_size = history_size;
		}

		void set_projection_lp(_double_t projection_lp) {
			this->projection_lp = projection_lp;
		}

		void set_gamma_factor(_double_t gamma_factor) {
			this->gamma_factor = gamma_factor;
		}

		void set_pcg_tol(_double_t pcg_tol) {
			this->pcg_tol = pcg_tol;
		}

		void set_rho_upper_limit(_double_t rho_upper_limit) {
			this->rho_upper_limit = rho_upper_limit;
		}

		void set_rel_tol(_double_t rel_tol) {
			this->rel_tol = rel_tol;
		}

		void set_pcg_maxiters(int pcg_maxiters) {
			this->pcg_maxiters = pcg_maxiters;
		}

		/* Setting the preconditioner for the conjugate gradient */
		void set_preconditioner(PreconditionerType preconditioner) {
			this->preconditioner = preconditioner;
		}

		void set_log_file(std::string st) {
			does_log = 1;
			this->log_file_path = st;
		}

};


#endif //CPP_LPBOXADMMSOLVER_H
