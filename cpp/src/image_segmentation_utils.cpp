#include "image_segmentation_utils.h"
#include <limits>


/* input = [a_1 a_2 ... a_n], output = [a_1; a_2; ... a_n] where a_i are column vectors */
void vectorize(const DenseMatrix &input, DenseVector &output) {
	int numRows = input.rows();
	int numCols = input.cols();
	output = DenseVector(numRows * numCols);
	for (int i = 0; i < numCols; i++) {
		output.block(i * numRows, 0, numRows, 1).array() = input.block(0, i, numRows, 1).array();
	}
}

void get_unary_cost(const DenseMatrix &_nodes, _double_t sigma, _double_t b, _double_t f1, _double_t f2, DenseMatrix &unary_cost) {
	DenseVector nodes;
	vectorize(_nodes, nodes);
	const _double_t c = std::log(2.0 * M_PI) / 2.0 + std::log(sigma);
	DenseVector alpha_b = (nodes.array() - b).pow(2.0) / (2 * sigma * sigma) + c;
	DenseVector aa = (-(nodes.array() - f1).pow(2.0) / (2 * sigma * sigma)).exp() 
		+ (-(nodes.array() - f2).pow(2) / (2 * sigma * sigma)).exp();
	DenseVector alpha_f = - (aa.array() + std::numeric_limits<_double_t>::epsilon()).log() + c + std::log(2.0);
	
	unary_cost = DenseMatrix(2, alpha_b.rows());
	unary_cost.block(0, 0, 1, alpha_b.rows()) = alpha_b.transpose();
	unary_cost.block(1, 0, 1, alpha_f.rows()) = alpha_f.transpose();
}

DenseMatrix get_row(const DenseMatrix &mat, int row_index) {
	int num_cols = mat.cols();
	DenseMatrix res = mat.block(row_index, 0, 1, num_cols);
	return res;
}

void repmat(const DenseMatrix &mat, int numRows, int numCols, DenseMatrix &result) {

	result = DenseMatrix(numRows * mat.rows(), numCols * mat.cols());
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			result.block(i * mat.rows(), j * mat.cols(), mat.rows(), mat.cols()) = mat;
		}
	}
}

/* the output is [A; B] */
void verticle_concat(const DenseMatrix &A, const DenseMatrix &B, DenseMatrix &output) {
	output = DenseMatrix(A.rows() + B.rows(), A.cols());
	output.block(0, 0, A.rows(), A.cols()) = A;
	output.block(A.rows(), 0, B.rows(), B.cols()) = B;
}

void verticle_concat(const DenseVector &A, const DenseVector &B, DenseVector &output) {
	output = DenseVector(A.rows() + B.rows());
	output.block(0, 0, A.rows(), 1) = A;
	output.block(A.rows(), 0, B.rows(), 1) = B;
}
	
void get_offsets(int k, DenseMatrix &offsets) {
	int edgeLen = 2 * k + 1;
	offsets = DenseMatrix(edgeLen * edgeLen, 2);
	DenseMatrix row_offsets = DenseMatrix(edgeLen, edgeLen);
	for (int i = 0; i < edgeLen; i++) {
		row_offsets.block(i, 0, 1, edgeLen).array() = i - k;
	}
	DenseMatrix col_offsets = row_offsets.transpose();

	for (int i = 0; i < edgeLen; i++) {
		offsets.block(i * edgeLen, 0, edgeLen, 1) = row_offsets.block(0, i, edgeLen, 1);
		offsets.block(i * edgeLen, 1, edgeLen, 1) = col_offsets.block(0, i, edgeLen, 1);
	}
}

void meshgrid(int x, int y, DenseMatrix &X, DenseMatrix &Y) {
	X = DenseMatrix(y, x);
	Y = DenseMatrix(y, x);
	for (int i = 0; i < x; i++) {
		X.block(0, i, y, 1).array() = i + 1;
	}

	for (int i = 0; i < y; i++) {
		Y.block(i, 0, 1, x).array() = i + 1;
	}
}

struct pair {
	int idx1;
	int idx2;
};

void generate_pixel_pairs(int nrows, int ncols, int k, DenseIntMatrix &res) {

	//printf("matrix rows: %d", ((2 * k + 1) * (2 * k + 1) - 1) * nrows * ncols);
	DenseIntMatrix tmp(((2 * k + 1) * (2 * k + 1) - 1) * nrows * ncols, 2);
	//std::vector<pair> vec;
	//vec.reserve();
	int counter = 0;
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			for (int a = -k; a < k + 1; a++) {
				for (int b = -k; b < k + 1; b++) {
					if (a != b && i + a >= 0 && i + a < nrows && j + b >= 0 && j + b < ncols) {
						//vec.push_back({i * ncols + j + 1 , (i + a) * ncols + (j + b + 1)});
						tmp(counter, 0) = i * ncols + j + 1;
						tmp(counter, 1) = (i + a) * ncols + (j + b + 1);
						counter += 1;
					}
				}
			}
		}
	}
	res = tmp.block(0, 0, counter, 2).matrix();
	
	//for (int i = 0; i < vec.size(); i++) {
		//res(i, 0) = vec[i].idx1;
		//res(i, 1) = vec[i].idx2;
	//}
}

void get_binary_cost(const DenseMatrix &image, SparseMatrix &binary_cost_matrix) {

	int nrows = image.rows();
	int ncols = image.cols();
	int num_pixels = nrows * ncols;

	DenseVector image_vec;
	vectorize(image, image_vec);
	_double_t sigma = std::sqrt((image_vec.array() - image_vec.mean()).square().sum()/(image_vec.size()-1));;

    uint8_t *indicator = new uint8_t[num_pixels];
    memset(indicator, 1, num_pixels);

	DenseIntMatrix pairs;
	generate_pixel_pairs(nrows, ncols, 1, pairs);
	DenseVector I1(pairs.rows());
	DenseVector I2(pairs.rows());

	for (int i = 0; i < pairs.rows(); i++) {
		I1(i) = image((pairs(i, 0) - 1) % nrows, (pairs(i, 0) - 1) / nrows);
		I2(i) = image((pairs(i, 1) - 1) % nrows, (pairs(i, 1) - 1) / nrows);
	}
	
	DenseVector diff_vec = (I1.array() - I2.array()).pow(2.0) / sigma;
	
	diff_vec.array() = (-diff_vec.array()).exp();

	binary_cost_matrix = SparseMatrix(num_pixels, num_pixels);
	
	std::vector<Triplet> triplets;
	for (int i = 0; i < pairs.rows(); i++) {
		int row = pairs(i, 0);
		int col = pairs(i, 1);
		triplets.push_back({row - 1, col - 1, std::round(3 * diff_vec(i))});
		if (row == col) {
			//printf("row: %d\n", row);
			indicator[row - 1] = 0;
		}
	}

	for (int i = 0; i < num_pixels; i++) {
		if (indicator[i] == 1) {
			//printf("row: %d\n", i);
			//printf("Encountering zero diagonal element\n");
			triplets.push_back({i, i, 0.0});
		}
	}

	binary_cost_matrix.setFromTriplets(triplets.begin(), triplets.end());
	delete[] indicator;

}

void get_A_b_from_cost(const DenseMatrix &unary_cost, const SparseMatrix &binary_cost, SparseMatrix &A, DenseVector &b) {
	int n = binary_cost.rows();
	
	DenseMatrix U1 = get_row(unary_cost, 0);
	DenseMatrix U2 = get_row(unary_cost, 1);

	b = (U2 - U1).eval().transpose();
	
	A = - binary_cost;
	DenseVector ones = DenseVector::Ones(binary_cost.cols());
	DenseVector We = - A * ones;
	A.diagonal().array() += We.array();
	A = 2 * A;
}

