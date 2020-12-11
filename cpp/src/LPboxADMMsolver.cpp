//
// Created by squintingsmile on 2020/9/22.
//

#include <cstring>
#include <cmath>
#include <chrono>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include "LPboxADMMsolver.h"

using namespace Eigen::internal;
using std::abs;
using std::sqrt;
Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
		", ", ", ", "", "", " << ", ";");

/* The original function is imported from the Eigen library. Added clocks to monitor times */
void _conjugate_gradient(const SparseMatrix& mat, const DenseVector& rhs, DenseVector& x,
		const Eigen::DiagonalPreconditioner<_double_t>& precond, int& iters, typename DenseVector::RealScalar& tol_error)
{

	typedef typename DenseVector::RealScalar RealScalar;
	typedef typename DenseVector::Scalar Scalar;
	typedef Eigen::Matrix<Scalar,Dynamic,1> VectorType;

	RealScalar tol = tol_error;
	int maxIters = iters;

	int n = mat.cols();

	VectorType residual = rhs - mat * x; //initial residual

	RealScalar rhsNorm2 = rhs.squaredNorm();

	if(rhsNorm2 == 0)
	{
		x.setZero();
		iters = 0;
		tol_error = 0;
		return;
	}

	const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
	RealScalar threshold = Eigen::numext::maxi(RealScalar(tol*tol*rhsNorm2),considerAsZero);
	RealScalar residualNorm2 = residual.squaredNorm();

	if (residualNorm2 < threshold)
	{
		iters = 0;
		tol_error = sqrt(residualNorm2 / rhsNorm2);
		return;
	}

	VectorType p(n);
	p = precond.solve(residual);      // initial search direction

	VectorType z(n), tmp(n);
	RealScalar absNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
	int i = 0;
	while(i < maxIters)
	{
		tmp.noalias() = mat * p;                    // the bottleneck of the algorithm

		Scalar alpha = absNew / p.dot(tmp);         // the amount we travel on dir

		x += alpha * p;                             // update solution

		residual -= alpha * tmp;                    // update residual

		residualNorm2 = residual.squaredNorm();

		if(residualNorm2 < threshold) {
			i++;
			break;
		}

		z = precond.solve(residual);                // approximately solve for "A z = residual"

		RealScalar absOld = absNew;
		absNew = Eigen::numext::real(residual.dot(z));     // update the absolute value of r
		RealScalar beta = absNew / absOld;          // calculate the Gram-Schmidt value used to create the new search direction
		p = z + beta * p;                           // update search direction
		i++;
	}
	tol_error = sqrt(residualNorm2 / rhsNorm2);
	iters = i;
	return;
}

/* The .noalias() function explicity states that result of the matrix multiplication can be written
 * on lhs directly without allocating an empty vector first and store the result, then assigning
 * it to lhs. This hopefully will save some time
 */
void mat_mul_vec(const SparseMatrix &mat, const DenseVector &vec, DenseVector& res) {
	if (mat.cols() == mat.rows() || &vec != &res) {
		res.noalias() = mat * vec;
	} else {
		res = mat * vec;
	}
}

/* This is used to avoid direct matrix multiplication. Suppose that we want to calculate A^TAx where A is large
 * so calculating A^TAx is intractable. Then we calculate Ax then calculate A^T(Ax) to reduce the number of operations.
 * In this vector mat_expression, each entry is a vector of sparse matrices. For each entry, the result is calculated by
 * entry[n] ... entry[2] * (entry[1] * (entry[0] * x)) and the results calculated by different entries are sum together
 */
void calculate_mat_expr_multiplication(const std::vector<std::vector<const SparseMatrix*>> &mat_expressions,
		DenseVector &x, DenseVector &result, DenseVector &temp_vec) {

	if (mat_expressions.size() == 0) {
		printf("Empty matrix expression when doing multiplication\n");
		return;
	}

	DenseVector *original_x;
	if (&x == &result && mat_expressions.size() > 1) {
		original_x = new DenseVector(x);
	} else {
		original_x = &x;
	}

	{
		const std::vector<const SparseMatrix*> &expression = mat_expressions[0];
		int counter = 0;
		for (const SparseMatrix *mat : expression) {
			if (counter == 0) {
				mat_mul_vec(*mat, x, result);
			} else {
				mat_mul_vec(*mat, result, result);
			}
			counter += 1;
		}
	}

	if (mat_expressions.size() > 1) {
		for (auto iter = std::next(mat_expressions.begin(), 1); iter != mat_expressions.end(); iter++) {
			int counter = 0;
			for (const SparseMatrix *mat : *iter) {
				if (counter == 0) {
					mat_mul_vec(*mat, *original_x, temp_vec);
				} else {
					mat_mul_vec(*mat, temp_vec, temp_vec);
				}
				counter += 1;
			}

			result += temp_vec;
		}

		if (&x == &result) {
			delete original_x;
		}
	}
}


/* The original function is imported from the Eigen library. Passing the two temporary vectors to save the time of
 * allocating vectors
 */

void _conjugate_gradient(const std::vector<std::vector<const SparseMatrix*>> &mat_expressions, const DenseVector& rhs,
		DenseVector& x, const Eigen::DiagonalPreconditioner<_double_t>& precond, int& iters,
		typename DenseVector::RealScalar& tol_error, DenseVector &temp_vec_for_cg, DenseVector &temp_vec_for_mat_mul) {

	typedef typename DenseVector::RealScalar RealScalar;
	typedef typename DenseVector::Scalar Scalar;
	typedef Eigen::Matrix<Scalar,Dynamic,1> VectorType;

	RealScalar tol = tol_error;
	int maxIters = iters;

	int n = mat_expressions[0][0]->cols();


	calculate_mat_expr_multiplication(mat_expressions, x, temp_vec_for_cg, temp_vec_for_mat_mul);
	VectorType residual = rhs - temp_vec_for_cg; //initial residual

	RealScalar rhsNorm2 = rhs.squaredNorm();

	if(rhsNorm2 == 0) {
		x.setZero();
		iters = 0;
		tol_error = 0;
		return;
	}

	const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
	RealScalar threshold = Eigen::numext::maxi(RealScalar(tol*tol*rhsNorm2),considerAsZero);
	RealScalar residualNorm2 = residual.squaredNorm();

	if (residualNorm2 < threshold)
	{
		iters = 0;
		tol_error = sqrt(residualNorm2 / rhsNorm2);
		return;
	}
	VectorType p(n);
	p = precond.solve(residual);      // initial search direction

	VectorType z(n), tmp(n);
	RealScalar absNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
	int i = 0;
	while(i < maxIters)
	{
		calculate_mat_expr_multiplication(mat_expressions, p, tmp, temp_vec_for_mat_mul);

		Scalar alpha = absNew / p.dot(tmp);         // the amount we travel on dir

		x += alpha * p;                             // update solution

		residual -= alpha * tmp;                    // update residual

		residualNorm2 = residual.squaredNorm();

		if(residualNorm2 < threshold) {
			i++;
			break;
		}

		z = precond.solve(residual);                // approximately solve for "A z = residual"

		RealScalar absOld = absNew;
		absNew = Eigen::numext::real(residual.dot(z));     // update the absolute value of r
		RealScalar beta = absNew / absOld;          // calculate the Gram-Schmidt value used to create the new search direction
		p = z + beta * p;                           // update search direction
		i++;
	}
	tol_error = sqrt(residualNorm2 / rhsNorm2);
	iters = i;
	return;
}


_double_t LPboxADMMsolver::std_dev(std::vector<_double_t>& arr, size_t begin, size_t end) {
	_double_t mean = 0;
	_double_t std_deviation = 0;
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

void LPboxADMMsolver::project_vec_greater_than(DenseVector &x, DenseVector &res, const int &greater_val, const int &set_val) {
	int len = x.size();
	for (int i = 0; i < len; i++) {
		res[i] = x[i] > greater_val ? set_val : x[i];
	}
}

void LPboxADMMsolver::project_vec_less_than(DenseVector &x, DenseVector &res, const int &less_val, const int &set_val) {
	int len = x.size();
	for (int i = 0; i < len; i++) {
		res(i) = x(i) < less_val ? set_val : x(i);
	}
}

void LPboxADMMsolver::project_box(int n, const double *x, double *y) {
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

/* TODO: Utilize eigen .array() functions to speed up this operation */
void LPboxADMMsolver::project_box(int n, const DenseVector &x, DenseVector &y) {
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

void LPboxADMMsolver::project_shifted_Lp_ball(int n, const DenseVector &x, int p, DenseVector &y) {
	y.array() = x.array() - 0.5;
	_double_t normp_shift = y.norm();

	y.array() = y.array() * std::pow(n, 1.0 / p) / (2 * normp_shift) + 0.5;
}

_double_t LPboxADMMsolver::compute_cost(const DenseVector& x, const SparseMatrix& A, 
		const DenseVector& b, DenseVector &temp_vec_for_mat_mul) {
	mat_mul_vec(A, x, temp_vec_for_mat_mul);
	double val = x.transpose().dot(temp_vec_for_mat_mul);
	double val2 = b.dot(x);
	return val + val2;
}

_double_t LPboxADMMsolver::compute_cost(const DenseVector& x, const SparseMatrix& A, const DenseVector& b) {
	double val = x.transpose().dot(A * x);
	double val2 = b.dot(x);
	return val + val2;
}

_double_t LPboxADMMsolver::compute_std_obj(std::vector<_double_t> obj_list, int history_size) {
	size_t s = obj_list.size();
	_double_t std_obj;
	if (s <= history_size) {
		std_obj = std_dev(obj_list, 0, s);
	} else {
		std_obj = std_dev(obj_list, s - history_size, s);
	}

	return std_obj / std::abs(obj_list[s-1]);
}

void LPboxADMMsolver::ADMM_bqp_unconstrained_init() {
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

void LPboxADMMsolver::ADMM_bqp_linear_eq_init() {
	stop_threshold = 1e-4;
	std_threshold = 1e-6;
	gamma_val = 1.6;
	gamma_factor = 0.95;
	rho_change_step = 5;
	max_iters = 5e3;
	initial_rho = 1;
	history_size = 3;
	learning_fact = 1 + 5.0 / 100;
	pcg_tol = 1e-4;
	pcg_maxiters = 1e3;
	rel_tol = 5e-5;
	projection_lp=2;
}

void LPboxADMMsolver::ADMM_bqp_linear_ineq_init() {
	stop_threshold = 1e-4;
	std_threshold = 1e-6;
	gamma_val = 1.6;
	gamma_factor = 0.95;
	rho_change_step = 5;
	max_iters = 1e4;
	initial_rho = 25;
	history_size = 3;
	learning_fact = 1 + 1.0 / 100;
	pcg_tol = 1e-4;
	pcg_maxiters = 1e3;
	rel_tol = 5e-5;
	projection_lp = 2;

}

void LPboxADMMsolver::ADMM_bqp_linear_eq_and_uneq_init() {
	stop_threshold = 1e-4;
	std_threshold = 1e-6;
	gamma_val = 1.6;
	gamma_factor = 0.95;
	rho_change_step = 5;
	max_iters = 1e4;
	initial_rho = 25;
	history_size = 3;
	learning_fact = 1 + 1.0 / 100;
	pcg_tol = 1e-4;
	pcg_maxiters = 1e3;
	rel_tol = 5e-5;
	projection_lp = 2;
}


int LPboxADMMsolver::ADMM_bqp(MatrixInfo matrix_info, SolverInstruction instruction, 
		Solution& sol) {


	const SparseMatrix *A_ptr, *C_ptr, *E_ptr;
	const DenseVector *b_ptr, *d_ptr, *f_ptr;

	A_ptr = matrix_info.A;
	b_ptr = matrix_info.b;

	SparseMatrix C_transpose, E_transpose, _2A_plus_rho1_rho2, rho3_C_transpose, rho4_E_transpose;

	DenseVector y3, z3, z4, Csq_diag, Esq_diag;

	int l = matrix_info.l;
	int m = matrix_info.m;
	int n = matrix_info.n;

	FILE *fp;
	if (does_log) {
		fp = fopen(log_file_path.c_str(), "w+");
	}

	DenseVector x_sol    = DenseVector(n);
	DenseVector y1       = DenseVector(n);
	DenseVector y2       = DenseVector(n);
	DenseVector z1       = DenseVector(n);
	DenseVector z2       = DenseVector(n);
	DenseVector prev_idx = DenseVector(n);
	DenseVector best_sol = DenseVector(n);
	DenseVector temp_vec = DenseVector(n);
	DenseVector cur_idx  = DenseVector(n);
	SparseMatrix temp_mat = SparseMatrix(n, n);

	/* The temp_vec_for_cg and temp_vec_for_mat_mat are vector that used to store the temporary vectors
	 * in the algorithm to prevent allocating new vectos. Adding this shall has around 5% of performance
	 * improvement (although might reduce the performance under some tasks by 10% since it increases the cost 
	 * of other operations besides matrix multiplication, such as vector addition and vector multiplication
	 * for reasons unknown)
	 */
	DenseVector temp_vec_for_cg = DenseVector(n);
	DenseVector temp_vec_for_mat_mul = DenseVector(n);
	_double_t cur_obj;
	bool rhoUpdated = true;

	x_sol = *(matrix_info.x0);

	z1.array() = 0;
	z2.array() = 0;
	cur_idx.array() = 0;


	Eigen::DiagonalPreconditioner<_double_t> diagonalPreconditioner;

	_double_t rho1 = initial_rho;
	_double_t rho2 = initial_rho;
	_double_t rho3 = initial_rho;
	_double_t rho4 = initial_rho;
	_double_t prev_rho1 = rho1;
	_double_t prev_rho2 = rho2;
	_double_t prev_rho3 = rho3;
	_double_t prev_rho4 = rho4;

	std::vector<_double_t> obj_list; /* Stores the objective value calculated during each iteration */
	_double_t std_obj = 1;

	_double_t cvg_test1;
	_double_t cvg_test2;
	_double_t rho_change_ratio;

	if (instruction.update_y3) {
		y3 = DenseVector(l);
	}

	if (instruction.update_z3) {
		z3 = DenseVector(m);
		z3.array() = 0;
	}

	if (instruction.update_z4) {
		z4 = DenseVector(l);
		z4.array() = 0;
	}


	/* If this task contains equality constraints, initialize C and d and relevant matrices */
	if (instruction.problem_type & equality) {
		C_ptr = matrix_info.C;
		C_transpose = (*C_ptr).transpose();
		rho3_C_transpose = rho3 * C_transpose;
		d_ptr = matrix_info.d;
	}

	/* If this task contains inequality constraints, initialize E and f and relevant matrices */
	if (instruction.problem_type & inequality) {
		E_ptr = matrix_info.E;
		E_transpose = (*E_ptr).transpose();
		rho4_E_transpose = rho4 * E_transpose;
		f_ptr = matrix_info.f;
	}

	/* Storing the matrix 2 * A + (rho1 + rho2) * I to save calculation time */
	_2A_plus_rho1_rho2 = 2 * (*A_ptr);
	_2A_plus_rho1_rho2.diagonal().array() += rho1 + rho2;

	SparseMatrix preconditioner_diag_mat(n, n); /* The diagonal matrix used by the preconditioner */
	preconditioner_diag_mat.reserve(n);
	std::vector<Triplet> preconditioner_diag_mat_triplets;
	preconditioner_diag_mat_triplets.reserve(n);
	/* The diagonal elements in the original expression evaluated by the preconditioner is given by the
	 * diagonal elements of 2 * _A + (rho1 + rho2) * I + rho3 * _C^T * _C
	 */
	for (int i = 0; i < n; i++) {
		preconditioner_diag_mat_triplets.push_back(Triplet(i, i, 0));
	}
	preconditioner_diag_mat.setFromTriplets(preconditioner_diag_mat_triplets.begin(), preconditioner_diag_mat_triplets.end());
	preconditioner_diag_mat.diagonal().array() = _2A_plus_rho1_rho2.diagonal().array();
	preconditioner_diag_mat.makeCompressed();


	/* The matrix expression for the unconstraint case is 2 * _A + (rho1 + rho2) * I */
	std::vector<std::vector<const SparseMatrix*>> matrix_expressions;
	matrix_expressions.emplace_back();
	matrix_expressions.back().push_back(&_2A_plus_rho1_rho2);

	/* If the problem has equality constraint, add rho3 * C^T * C to the matrix expression */
	if (instruction.problem_type & equality) {
		matrix_expressions.emplace_back();
		matrix_expressions.back().push_back(C_ptr);
		matrix_expressions.back().push_back(&rho3_C_transpose);

		/* Calculating the diagonal elements of Csq for preconditioner */
		Csq_diag = DenseVector(n); 		
		Csq_diag.setZero();
		for(int j = 0; j < C_transpose.outerSize(); ++j) {
			typename SparseMatrix::InnerIterator it(C_transpose, j);
			while (it) {
				if(it.value() != 0.0) {
					Csq_diag[j] += it.value() * it.value();
				}
				++it;
			}
		}
		preconditioner_diag_mat.diagonal().array() += rho3 * Csq_diag.array();
	}

	/* If the problem has equality constraint, add rho4 * E^T * E to the matrix expression */
	if (instruction.problem_type &inequality) {
		matrix_expressions.emplace_back();
		matrix_expressions.back().push_back(E_ptr);
		matrix_expressions.back().push_back(&rho4_E_transpose);

		/* Calculating the diagonal elements of Esq for preconditioner */
		Esq_diag = DenseVector(n);
		Esq_diag.setZero();
		for(int j = 0; j < E_transpose.outerSize(); ++j) {
			typename SparseMatrix::InnerIterator it(E_transpose, j);
			while (it) {
				if(it.value() != 0.0) {
					Esq_diag[j] += it.value() * it.value();
				}
				++it;
			}
		}
		preconditioner_diag_mat.diagonal().array() += rho4 * Esq_diag.array();
	}

	y1 = x_sol;
	y2 = x_sol;
	if (instruction.update_y3) {
		y3 = *f_ptr - *E_ptr * x_sol;
	}

	prev_idx = (x_sol.array() >= 0.5).matrix().cast<_double_t>();
	best_sol = x_sol;

	_double_t best_bin_obj = compute_cost(x_sol, *A_ptr, *b_ptr, temp_vec_for_mat_mul);

	if (does_log) {
		fprintf(fp, "Initial state\n");
        fprintf(fp, "norm of x_sol: %lf\n", x_sol.norm());
        fprintf(fp, "norm of b: %lf\n", (*b_ptr).norm());
        fprintf(fp, "norm of y1: %lf\n", y1.norm());
        fprintf(fp, "norm of y2: %lf\n", y2.norm());
        if (instruction.update_y3) {
            fprintf(fp, "norm of y3: %lf\n", y3.norm());
        }

        fprintf(fp, "norm of z1: %lf\n", z1.norm());
        fprintf(fp, "norm of z2: %lf\n", z2.norm());

        if (instruction.update_z3) {
            fprintf(fp, "norm of z3: %lf\n", z3.norm());
        }

        if (instruction.update_z4) {
            fprintf(fp, "norm of z4: %lf\n", z4.norm());
        }

        fprintf(fp, "norm of cur_idx: %lf\n", cur_idx.norm());
        fprintf(fp, "rho1: %lf\n", rho1);
        fprintf(fp, "rho2: %lf\n", rho2);
        fprintf(fp, "rho3: %lf\n", rho3);
        fprintf(fp, "rho4: %lf\n", rho4);
        fprintf(fp, "-------------------------------------------------\n");
	}

	std::chrono::steady_clock::time_point start, end;
	start = std::chrono::steady_clock::now();


	long time_elapsed = 0;
	for (int iter = 0; iter < max_iters; iter++) {
		if (does_log) {
			fprintf(fp, "Iteration: %d\n", iter);
		}
		temp_vec = x_sol + z1 / rho1;

		/* Project vector on [0, 1] box */
		project_box(n, temp_vec, y1);

		temp_vec = x_sol + z2 / rho2;

		/* Project vector on shifted lp box */
		project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

		if (instruction.update_y3) {
			mat_mul_vec(*E_ptr, x_sol, temp_vec_for_mat_mul);
			y3 = *f_ptr - temp_vec_for_mat_mul - z4 / rho4;
			project_vec_less_than(y3, y3, 0, 0); 
		}

		/* If the iteration is nonzero and it divides rho_change_step, it means
		 * that the rho updated in the last iteration
		 */
		if (iter != 0 && rhoUpdated) {
			/* Note that we need the previous rho to update in order to get the difference between
			 * the updated matrix and the not updated matrix. Another possible calculation is by
			 * calculating rho * rho_change_ration / learning_fact
			 */
			_2A_plus_rho1_rho2.diagonal().array() += rho_change_ratio * (prev_rho1 + prev_rho2);
			if (instruction.problem_type != unconstrained) {
				preconditioner_diag_mat.diagonal().array() += rho_change_ratio * (prev_rho1 + prev_rho2);
			}

			if (instruction.update_rho3) {
				preconditioner_diag_mat.diagonal().array() += rho_change_ratio * prev_rho3 * Csq_diag.array();
				rho3_C_transpose = learning_fact * rho3_C_transpose;	
			}

			if (instruction.update_rho4) {
				preconditioner_diag_mat.diagonal().array() += rho_change_ratio * prev_rho4 * Esq_diag.array();
				rho4_E_transpose = learning_fact * rho4_E_transpose;
			}
		}



		/* If the problem in unconstrained, the rhs vector is 
		 * rho1 * y1 + rho2 * y2 - (b + z1 + z2)
		 */
		if (instruction.problem_type == unconstrained) {
			temp_vec = rho1 * y1 + rho2 * y2 - (*b_ptr + z1 + z2);
		}

		/* If the problem in equality, the rhs vector is 
		 * rho1 * y1 + rho2 * y2 + rho3 * C^T * d - (b + z1 + z2 + C^T * z3)
		 */
		if (instruction.problem_type == equality) {
			temp_vec = rho1 * y1 + rho2 * y2 - (*b_ptr + z1 + z2);
			mat_mul_vec(rho3_C_transpose, *d_ptr, temp_vec_for_mat_mul);
			temp_vec += temp_vec_for_mat_mul;
			mat_mul_vec(C_transpose, z3, temp_vec_for_mat_mul);
			temp_vec -= temp_vec_for_mat_mul;
		}

		/* If the problem in equality, the rhs vector is 
		 * rho1 * y1 + rho2 * y2 + rho4 * E^T * (f - y3) - (b + z1 + z2 + E^T * z4)
		 */
		if (instruction.problem_type == inequality) {
			temp_vec = rho1 * y1 + rho2 * y2 - (*b_ptr + z1 + z2);
			mat_mul_vec(rho4_E_transpose, *f_ptr - y3, temp_vec_for_mat_mul);
			temp_vec += temp_vec_for_mat_mul;
			mat_mul_vec(E_transpose, z4, temp_vec_for_mat_mul);
			temp_vec -= temp_vec_for_mat_mul;
		}

		/* If the problem in equality, the rhs vector is 
		 * rho1 * y1 + rho2 * y2 + rho3 * C^T * d + rho4 * E^T * (f - y3) - (b + z1 + z2 + C^T * z3 + E^T * z4)
		 */
		if (instruction.problem_type == equality_and_inequality) {
			temp_vec = rho1 * y1 + rho2 * y2 - (*b_ptr + z1 + z2);
			mat_mul_vec(rho3_C_transpose, *d_ptr, temp_vec_for_mat_mul);
			temp_vec += temp_vec_for_mat_mul;
			mat_mul_vec(C_transpose, z3, temp_vec_for_mat_mul);
			temp_vec -= temp_vec_for_mat_mul;
			mat_mul_vec(rho4_E_transpose, *f_ptr - y3, temp_vec_for_mat_mul);
			temp_vec += temp_vec_for_mat_mul;
			mat_mul_vec(E_transpose, z4, temp_vec_for_mat_mul);
			temp_vec -= temp_vec_for_mat_mul;
		}


		/* Explicit version of conjugate gradient used for profiling */
		if (rhoUpdated) {
			if (instruction.problem_type != unconstrained) {
				diagonalPreconditioner.compute(preconditioner_diag_mat);
			} else {
				diagonalPreconditioner.compute(_2A_plus_rho1_rho2);
			}
			rhoUpdated = false;
		}
		_double_t tol = pcg_tol;
		x_sol = y1;
		int maxiter = pcg_maxiters;
		_conjugate_gradient(matrix_expressions, temp_vec, x_sol, diagonalPreconditioner, 
				maxiter, tol, temp_vec_for_cg, temp_vec_for_mat_mul);

		if (does_log) {
			fprintf(fp, "Conjugate gradient stops after %d iterations\n", maxiter);
			fprintf(fp, "Conjugate gradient stops with residual %lf\n", tol);
		}

		z1 = z1 + gamma_val * rho1 * (x_sol - y1);
		z2 = z2 + gamma_val * rho2 * (x_sol - y2);
		if (instruction.update_z3) {
			z3 = z3 + gamma_val * rho3 * (*C_ptr * x_sol - *d_ptr);
		}

		if (instruction.update_z4) {
			z4 = z4 + gamma_val * rho4 * ((*E_ptr) * x_sol + y3 - (*f_ptr));
		}


		_double_t temp0 = std::max(x_sol.norm(), _double_t(2.2204e-16));
		cvg_test1 = (x_sol - y1).norm() / temp0;
		cvg_test2 = (x_sol - y2).norm() / temp0;
		if (cvg_test1 <= stop_threshold && cvg_test2 <= stop_threshold) {
			printf("iter: %d, stop_threshold: %.6f\n", iter, std::max(cvg_test1, cvg_test2));
			if (does_log) {
				fprintf(fp, "iter: %d, stop_threshold: %.6f\n", iter, std::max(cvg_test1, cvg_test2));
			}
			break;
		}

		if ((iter+1) % rho_change_step == 0) {
			prev_rho1 = rho1;
			prev_rho2 = rho2;
			rho1 = learning_fact * rho1;
			rho2 = learning_fact * rho2;

			if (instruction.update_rho3) {
				prev_rho3 = rho3;
				rho3 = learning_fact * rho3;
			}

			if (instruction.update_rho4) {
				prev_rho4 = rho4;
				rho4 = learning_fact * rho4;
			}

			gamma_val = std::max(gamma_val * gamma_factor, _double_t(1.0));
			rhoUpdated = true;
			rho_change_ratio = learning_fact - 1.0;
		}

		_double_t obj_val = compute_cost(x_sol, *A_ptr, *b_ptr);
		obj_list.push_back(obj_val);
		if (obj_list.size() >= history_size) {
			std_obj = compute_std_obj(obj_list, history_size);
		}
		if (std_obj <= std_threshold) {
			if (does_log) {
				fprintf(fp, "iter: %d, std_threshold: %.6f\n", iter, std_obj);
			}
			printf("iter: %d, std_threshold: %.6f\n", iter, std_obj);
			break;
		}

		cur_idx = (x_sol.array() >= 0.5).matrix().cast<_double_t>();
		prev_idx = cur_idx;
		cur_obj = compute_cost(prev_idx, *A_ptr, *b_ptr);

		if (best_bin_obj >= cur_obj) {
			best_bin_obj = cur_obj;
			best_sol = x_sol;
		}

		if (does_log) {
			fprintf(fp, "current objective: %lf\n", obj_val);
			fprintf(fp, "current binary objective: %lf\n", cur_obj);

			if (instruction.problem_type == equality || instruction.problem_type == equality_and_inequality) {
				fprintf(fp, "equality constraint violation: %lf\n", (*matrix_info.C * cur_idx - *matrix_info.d).norm() / x_sol.rows());
			}

			if (instruction.problem_type == inequality || instruction.problem_type == equality_and_inequality) {
				DenseVector diff = *matrix_info.E * cur_idx - *matrix_info.f;
				project_vec_less_than(diff, diff, 0, 0);
				fprintf(fp, "inequality constraint violation: %lf\n", diff.norm() / x_sol.rows());
			}

			fprintf(fp, "norm of x_sol: %lf\n", x_sol.norm());
			fprintf(fp, "norm of y1: %lf\n", y1.norm());
			fprintf(fp, "norm of y2: %lf\n", y2.norm());
			if (instruction.update_y3) {
				fprintf(fp, "norm of y3: %lf\n", y3.norm());
			}

			fprintf(fp, "norm of z1: %lf\n", z1.norm());
			fprintf(fp, "norm of z2: %lf\n", z2.norm());

			if (instruction.update_z3) {
				fprintf(fp, "norm of z3: %lf\n", z3.norm());
			}

			if (instruction.update_z4) {
				fprintf(fp, "norm of z4: %lf\n", z4.norm());
			}

			fprintf(fp, "norm of cur_idx: %lf\n", cur_idx.norm());
			fprintf(fp, "rho1: %lf\n", rho1);
			fprintf(fp, "rho2: %lf\n", rho2);
			if (instruction.update_rho3) {
				fprintf(fp, "rho3: %lf\n", rho3);
			}
			if (instruction.update_rho4) {
				fprintf(fp, "rho4: %lf\n", rho4);
			}
			fprintf(fp, "-------------------------------------------------\n");
		}
	}

	sol.x_sol = new DenseVector(x_sol);
	sol.y1 = new DenseVector(y1);
	sol.y2 = new DenseVector(y2);
	sol.best_sol = new DenseVector(best_sol);

	end = std::chrono::steady_clock::now();
	time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time elapsed: " << time_elapsed << "us" << std::endl;
	if (does_log) {
		fprintf(fp, "Time elapsed: %ldus\n", time_elapsed);
		fclose(fp);
	}
	sol.time_elapsed = time_elapsed;

	return 1;
}

int LPboxADMMsolver::ADMM_bqp_unconstrained(int n, const SparseMatrix &_A, const DenseVector &_b, 
		const DenseVector &x0, Solution& sol) {

	SolverInstruction solver_instruction;
	MatrixInfo matrix_info;

	matrix_info.x0 = &x0;
	matrix_info.A = &_A;
	matrix_info.b = &_b;
	matrix_info.n = n;

	solver_instruction.problem_type = unconstrained;
	solver_instruction.update_y3 = 0;
	solver_instruction.update_z3 = 0;
	solver_instruction.update_z4 = 0;
	solver_instruction.update_rho3 = 0;
	solver_instruction.update_rho4 = 0;

	ADMM_bqp(matrix_info, solver_instruction, sol);
	return 1;
}



int LPboxADMMsolver::ADMM_bqp_unconstrained(int n, _double_t *A, _double_t *b, _double_t *x0, 
		Solution& sol) {
	/* Initialize the parameters */
	auto _A = SparseMatrix(n, n);
	auto _b = DenseVector(n);
	auto _x0 = DenseVector(n);

	/* Generating the sparse matrix for A */
	std::vector<Triplet> triplet_list;
	for (int i = 0; i < n * n; i++) {
		int row = i % n;
		int col = i / n;

		if (A[i] == 0.0 && !(row == col)) {
			continue;
		}
		triplet_list.push_back(Triplet(row, col, A[i]));
	}
	_A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	memcpy(_b.data(), b, n * sizeof(_double_t));
	memcpy(_x0.data(), x0, n * sizeof(_double_t));

	int ret = ADMM_bqp_unconstrained(n, _A, _b, _x0, sol);
	return ret;
}


int LPboxADMMsolver::ADMM_bqp_linear_eq(int n, const SparseMatrix &_A, const DenseVector &_b, 
		const DenseVector &x0, int m, const SparseMatrix &_C, const DenseVector &_d, Solution& sol) {

	SolverInstruction solver_instruction;
	MatrixInfo matrix_info;

	matrix_info.x0 = &x0;
	matrix_info.A = &_A;
	matrix_info.b = &_b;
	matrix_info.n = n;
	matrix_info.C = &_C;
	matrix_info.d = &_d; 
	matrix_info.m = m;

	solver_instruction.problem_type = equality;
	solver_instruction.update_y3 = 0;
	solver_instruction.update_z3 = 1;
	solver_instruction.update_z4 = 0;

	/* If this value is set to be true, then the admm update will update
     * rho3 at each iteration with rate learning_fact. The default setting
     * is false, following the strategy in the clustering task */
	solver_instruction.update_rho3 = 0;
	solver_instruction.update_rho4 = 0;

	ADMM_bqp(matrix_info, solver_instruction, sol);
	return 1;
}


int LPboxADMMsolver::ADMM_bqp_linear_eq(int n, _double_t *A, _double_t *b, _double_t *x0, int m, 
		_double_t *C, _double_t *d, Solution& sol) {
	/* Initialize the parameters */
	auto _A = SparseMatrix(n, n);
	auto _b = DenseVector(n);
	auto _x0 = DenseVector(n);
	auto _C = SparseMatrix(m, n);
	auto _d = DenseVector(m);

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

	memcpy(_x0.data(), x0, n * sizeof(_double_t));
	memcpy(_b.data(), b, n * sizeof(_double_t));
	memcpy(_d.data(), d, m * sizeof(_double_t));

	int ret = ADMM_bqp_linear_eq(n, _A, _b, _x0, m, _C, _d, sol);
	return ret;
}

int LPboxADMMsolver::ADMM_bqp_linear_ineq(int n, const SparseMatrix &_A, const DenseVector &_b,
		const DenseVector &x0, int l, const SparseMatrix &_E, const DenseVector &_f, Solution& sol) {

	SolverInstruction solver_instruction;
	MatrixInfo matrix_info;

	matrix_info.x0 = &x0;
	matrix_info.A = &_A;
	matrix_info.b = &_b;
	matrix_info.n = n;
	matrix_info.E = &_E;
	matrix_info.f = &_f; 
	matrix_info.l = l;

	solver_instruction.problem_type = inequality;
	solver_instruction.update_y3 = 1;
	solver_instruction.update_z3 = 0;
	solver_instruction.update_z4 = 1;
	solver_instruction.update_rho3 = 0;
	solver_instruction.update_rho4 = 1;

	ADMM_bqp(matrix_info, solver_instruction, sol);
	return 1;
}

int LPboxADMMsolver::ADMM_bqp_linear_ineq(int n, _double_t *A, _double_t *b, _double_t *x0, 
		int l, _double_t *E, _double_t *f, Solution& sol) {
	/* Initialize the parameters */
	auto _A       = SparseMatrix(n, n);
	auto _b       = DenseVector(n);
	auto _x0	  = DenseVector(n);
	auto _E       = SparseMatrix(l, n);
	auto _f       = DenseVector(l);

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
	for (int i = 0; i < l * n; i++) {
		if (E[i] == 0.0) {
			continue;
		}
		int row = i % n;
		int col = i / n;
		triplet_list.push_back(Triplet(row, col, E[i]));
	}
	_E.setFromTriplets(triplet_list.begin(), triplet_list.end());

	memcpy(_x0.data(), x0, n * sizeof(_double_t));
	memcpy(_b.data(), b, n * sizeof(_double_t));
	memcpy(_f.data(), f, l * sizeof(_double_t));

	int ret = ADMM_bqp_linear_ineq(n, _A, _b, _x0, l, _E, _f, sol);
	return ret;
}



int LPboxADMMsolver::ADMM_bqp_linear_eq_and_uneq(int n, const SparseMatrix &_A, 
		const DenseVector &_b, const DenseVector &x0, int m, const SparseMatrix &_C, const DenseVector &_d, 
		int l, const SparseMatrix &_E, const DenseVector &_f, Solution& sol) {

	SolverInstruction solver_instruction;
	MatrixInfo matrix_info;

	matrix_info.x0 = &x0;
	matrix_info.A = &_A;
	matrix_info.b = &_b;
	matrix_info.n = n;
	matrix_info.C = &_C;
	matrix_info.d = &_d;
	matrix_info.m = m;
	matrix_info.E = &_E;
	matrix_info.f = &_f; 
	matrix_info.l = l;

	solver_instruction.problem_type = equality_and_inequality;
	solver_instruction.update_y3 = 1;
	solver_instruction.update_z3 = 1;
	solver_instruction.update_z4 = 1;
	solver_instruction.update_rho3 = 1;
	solver_instruction.update_rho4 = 1;

	ADMM_bqp(matrix_info, solver_instruction, sol);
	return 1;

}


/* Currently assuming that the input matrix are column majored */
int LPboxADMMsolver::ADMM_bqp_linear_eq_and_uneq(int n, _double_t *A, _double_t *b, 
		_double_t *x0, int m, _double_t *C, _double_t *d, int l, _double_t *E, _double_t *f, Solution& sol) {
	/* Initialize the parameters */
	auto _A       = SparseMatrix(n, n);
	auto _b       = DenseVector(n);
	auto _x0      = DenseVector(n);
	auto _C       = SparseMatrix(m, n);
	auto _d       = DenseVector(m);
	auto _E       = SparseMatrix(l, n);
	auto _f       = DenseVector(l);

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
	for (int i = 0; i < l * n; i++) {
		if (C[i] == 0.0) {
			continue;
		}
		int row = i % n;
		int col = i / n;
		triplet_list.push_back(Triplet(row, col, C[i]));
	}
	_C.setFromTriplets(triplet_list.begin(), triplet_list.end());

	memcpy(_x0.data(), x0, n * sizeof(_double_t));
	memcpy(_b.data(), b, n * sizeof(_double_t));
	memcpy(_d.data(), d, m * sizeof(_double_t));
	memcpy(_f.data(), f, l * sizeof(_double_t));

	int ret = ADMM_bqp_linear_eq_and_uneq(n, _A, _b, _x0, m, _C, _d, l, _E, _f, sol);
	return ret;
}

int LPboxADMMsolver::ADMM_bqp_unconstrained_legacy(int n, const SparseMatrix &_A, const DenseVector &_b,
		const DenseVector &x0, Solution& sol) {

	auto x_sol    = DenseVector(n);
	auto y1       = DenseVector(n);
	auto y2       = DenseVector(n);
	auto z1       = DenseVector(n);
	auto z2       = DenseVector(n);
	auto prev_idx = DenseVector(n);
	auto best_sol = DenseVector(n);
	auto temp_vec = DenseVector(n);
	auto temp_mat = SparseMatrix(n, n);
	auto cur_idx  = DenseVector(n);
	_double_t cur_obj;
	bool rhoUpdated = true;

	x_sol = x0;

	/* Initializing preconditioner for conjugate gradient */
	Eigen::DiagonalPreconditioner<_double_t> diagonalPreconditioner;

	_double_t rho1 = initial_rho;
	_double_t rho2 = initial_rho;
	_double_t prev_rho1 = rho1;
	_double_t prev_rho2 = rho2;
	std::vector<_double_t> obj_list; /* Stores the objective value calculated during each iteration */
	_double_t std_obj = 1;

	/* temp_mat stores the matrix that is used in the conjugate gradient step */
	temp_mat = 2 * _A;
	temp_mat.diagonal().array() += rho1 + rho2;
	temp_mat.makeCompressed();

	_double_t cvg_test1;
	_double_t cvg_test2;
	_double_t rho_change_ratio;

	y1 = x_sol;
	y2 = x_sol;

	FILE *fp;
	if (does_log) {
		fp = fopen(log_file_path.c_str(), "w+");
	}

	prev_idx = (x_sol.array() >= 0.5).matrix().cast<_double_t>();
	best_sol = x_sol;

	_double_t best_bin_obj = compute_cost(x_sol, _A, _b);

	if (does_log) {
		fprintf(fp, "Initial state\n");
        fprintf(fp, "norm of x_sol: %lf\n", x_sol.norm());
        fprintf(fp, "norm of y1: %lf\n", y1.norm());
        fprintf(fp, "norm of y2: %lf\n", y2.norm());

        fprintf(fp, "norm of z1: %lf\n", z1.norm());
        fprintf(fp, "norm of z2: %lf\n", z2.norm());


        fprintf(fp, "norm of cur_idx: %lf\n", cur_idx.norm());
        fprintf(fp, "rho1: %lf\n", rho1);
        fprintf(fp, "rho2: %lf\n", rho2);
        fprintf(fp, "-------------------------------------------------\n");
	}

	std::chrono::steady_clock::time_point start, end;
	start = std::chrono::steady_clock::now();


	long time_elapsed = 0;
	for (int iter = 0; iter < max_iters; iter++) {
		temp_vec = x_sol + z1 / rho1;

		if (does_log) {
			fprintf(fp, "Iteration: %d\n", iter);
		}
		/* Project vector on [0, 1] box to calculate y1 */
		project_box(n, temp_vec, y1);

		temp_vec = x_sol + z2 / rho2;

		/* Project vector on shifted lp box to calculate y2 */
		project_shifted_Lp_ball(n, temp_vec, projection_lp, y2);

		/* If the iteration is nonzero and it divides rho_change_step, it means
		 * that the rho updated in the last iteration and the matrix used by conjugate
		 * gradient should be updated
		 */
		if (iter != 0 && rhoUpdated) {
			temp_mat.diagonal().array() += (prev_rho1 + prev_rho2) * rho_change_ratio;
			temp_mat.makeCompressed();
		}

		/* Calculate the vector b in the conjugate gradient algorithm */
		temp_vec = rho1 * y1 + rho2 * y2 - (_b + z1 + z2);

		/* Explicit version of conjugate gradient used for profiling */

		/* Since the matrix used by the conjugate gradient changes only when after rho is updated
		 * we only need to recalculate the preconditioner if rho is updated in the last iteration
		 */
		if (rhoUpdated) {
			diagonalPreconditioner.compute(temp_mat);
			rhoUpdated = false;
		}

		x_sol = y1;
		_double_t tol = pcg_tol;
		int maxiter = pcg_maxiters;
		_conjugate_gradient(temp_mat, temp_vec, x_sol, diagonalPreconditioner, maxiter, tol);
		if (does_log) {
			fprintf(fp, "Conjugate gradient stops after %d iterations\n", maxiter);
			fprintf(fp, "Conjugate gradient stops with residual %lf\n", tol);
		}

		z1 = z1 + gamma_val * rho1 * (x_sol - y1);
		z2 = z2 + gamma_val * rho2 * (x_sol - y2);


		/* Testing the conditions to see if the algorithm converges */
		_double_t temp0 = std::max(x_sol.norm(), _double_t(2.2204e-16));
		cvg_test1 = (x_sol-y1).norm() / temp0;
		cvg_test2 = (x_sol-y2).norm() / temp0;
		if (cvg_test1 <= stop_threshold && cvg_test2 <= stop_threshold) {
			printf("iter: %d, stop_threshold: %.6f\n", iter, std::max(cvg_test1, cvg_test2));
			if (does_log) {
				fprintf(fp, "iter: %d, stop_threshold: %.6f\n", iter, std::max(cvg_test1, cvg_test2));
			}
			break;
		}

		/* Update the rho value every rho_change_step */
		if ((iter+1) % rho_change_step == 0) {
			prev_rho1 = rho1;
			prev_rho2 = rho2;
			rho1 = learning_fact * rho1;
			rho2 = learning_fact * rho2;
			gamma_val = std::max(gamma_val * gamma_factor, _double_t(1.0));
			rhoUpdated = true;
			rho_change_ratio = learning_fact - 1.0;
		}

		/* Computer the relaxed cost function (x is not binary)*/
		_double_t obj_val = compute_cost(x_sol,_A,_b);
		obj_list.push_back(obj_val);
		if (obj_list.size() >= history_size) {
			std_obj = compute_std_obj(obj_list, history_size);
		}
		if (std_obj <= std_threshold) {
			printf("iter: %d, std_threshold: %.6f\n", iter, std_obj);
			if (does_log) {
				fprintf(fp, "iter: %d, std_threshold: %.6f\n", iter, std_obj);
			}
			break;
		}

		/* Calculating the actual cost */
		cur_idx = (x_sol.array() >= 0.5).matrix().cast<_double_t>(); /* The value type of the vector should be double
																	  therefore needs casting the boolean vector to double */
		prev_idx = cur_idx;
		cur_obj = compute_cost(prev_idx, _A, _b);

		/* Setting the best binary solution */
		if (best_bin_obj >= cur_obj) {
			best_bin_obj = cur_obj;
			best_sol = x_sol;
		}

		if (does_log) {
			fprintf(fp, "current objective: %lf\n", obj_val);
			fprintf(fp, "current binary objective: %lf\n", cur_obj);
			fprintf(fp, "norm of x_sol: %lf\n", x_sol.norm());
			fprintf(fp, "norm of binary x_sol: %lf\n", cur_idx.norm());
			fprintf(fp, "norm of y1: %lf\n", y1.norm());
			fprintf(fp, "norm of y2: %lf\n", y2.norm());

			fprintf(fp, "norm of z1: %lf\n", z1.norm());
			fprintf(fp, "norm of z2: %lf\n", z2.norm());

			fprintf(fp, "rho1: %lf\n", rho1);
			fprintf(fp, "rho2: %lf\n", rho2);
			fprintf(fp, "-------------------------------------------------\n");
		}

	}
	sol.x_sol = new DenseVector(x_sol);
	sol.y1 = new DenseVector(y1);
	sol.y2 = new DenseVector(y2);
	sol.best_sol = new DenseVector(best_sol);

	end = std::chrono::steady_clock::now();
	time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time elapsed: " << time_elapsed << "us" << std::endl;
	sol.time_elapsed = time_elapsed;

	if (does_log) {
		fprintf(fp, "Time elapsed: %ldus\n", time_elapsed);
		fclose(fp);
	}
	return 1;
}
