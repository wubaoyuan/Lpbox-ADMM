//
// Created by squintingsmile on 2020/10/5.
//

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <string>
#include <limits>
#include "read_matrix_utils.h"

struct Trip {
	int row;
	int col;
	_double_t val;
};

void readDenseMat(FILE *fp, DenseMatrix &mat, int rowNum, int colNum) {
	double val;
	mat = DenseMatrix(rowNum, colNum);

	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			if (fscanf(fp, "%lf,", &val) == 0) {
				printf("error when reading dense matrix\n");
				exit(-1);
			}
			mat(i, j) = val;
		}
	}
}

void readDenseMatAsSparse(FILE *fp, SparseMatrix &smat, int rowNum, int colNum) {
	double val;
	DenseMatrix mat(rowNum, colNum);

	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			if (fscanf(fp, "%lf,", &val) == 0) {
				printf("error when reading dense matrix as sparse matrix\n");
				exit(-1);
			}
			mat(i, j) = val;
		}
	}
	mat.diagonal().array() += std::numeric_limits<_double_t>::epsilon();

	smat = mat.sparseView();
}

void readDenseVec(FILE *fp, DenseVector &vec, int vecLen) {
	for (int i = 0; i < vecLen; i++) {
		if (fscanf(fp, "%lf,", &(vec.data()[i])) == 0) {
			printf("error when reading dense vector\n");
			exit(-1);
		}
	}
}


void readSparseMat(FILE *fp, SparseMatrix &mat) {
	std::vector<Triplet> triplets;
	int row, col;
	double val;
	int max_row = 0;
	int max_col = 0;
	while(true) {
		if (fscanf(fp, "%d,%d,%lf\n", &row, &col, &val) != 3) {
			break;
		}
		if (row > max_row) {
			max_row = row;
		}

		if (col > max_col) {
			max_col = col;
		}
		triplets.push_back(Triplet(row - 1, col - 1, val));
	}
	mat = SparseMatrix(max_row, max_col);
	mat.setFromTriplets(triplets.begin(), triplets.end());
}

/* Ensure that the diagonal elements exist in the sparse matrix */
void readSparseMat_diag(FILE *fp, SparseMatrix &mat) {
	std::vector<Triplet> triplets;
	int row, col;
	double val;
	std::vector<int> indicators;
	int max_row = 0;
	int max_col = 0;
	while(true) {
		if (fscanf(fp, "%d,%d,%lf\n", &row, &col, &val) != 3) {
			break;
		}
		if (row == col) {
			indicators.push_back(row);
		}
		if (row > max_row) {
			max_row = row;
		}

		if (col > max_col) {
			max_col = col;
		}
		triplets.push_back(Triplet(row - 1, col - 1, val));
	}

	uint8_t *indicator = new uint8_t[max_row];
	memset(indicator, 0, max_row);
	for (int i = 0; i < indicators.size(); i++) {
		indicator[indicators[i] - 1] = 1;	
	}

	for (int i = 0; i < indicators.size(); i++) {
		if (indicator[i] == 0) {
			triplets.push_back(Triplet(i, i, 0.0));
		}
	}

	delete[] indicator;
	mat = SparseMatrix(max_row, max_col);
	mat.setFromTriplets(triplets.begin(), triplets.end());
}

/* For the clustering mission, generate the matrices from the laplacian */
void generateMatrixFromLaplacian(const DenseMatrix &laplacian, SparseMatrix &A, DenseVector &b, 
		SparseMatrix &C, DenseVector &d, int numNodes, int numClasses) {
	DenseMatrix _A(numNodes * numClasses, numNodes * numClasses);
	_A.setZero();
	DenseMatrix _C(numNodes + numClasses, numNodes * numClasses);
	_C.setZero();

	DenseMatrix allOneRowVector = DenseMatrix(1, numNodes);
	allOneRowVector.array() = allOneRowVector.array() + 1.0;
	for (int i = 0; i < numClasses; i++) {
		_A.block(i * numNodes, i * numNodes, numNodes, numNodes) = laplacian;
		_C.block(i, i * numNodes, 1, numNodes) = allOneRowVector;
		_C.block(numClasses, i * numNodes, numNodes, numNodes) = DenseMatrix::Identity(numNodes, numNodes);
	}

	A = _A.sparseView();
	C = _C.sparseView();
	b = DenseVector::Zero(numNodes * numClasses);
	d = DenseVector::Zero(numNodes + numClasses);
	d.block(0, 0, numClasses, 1).array() = (double) numNodes / (double) numClasses;
	d.block(numClasses, 0, numNodes, 1).array() = 1;

}

void writeVectorToFile(const DenseVector &vec, FILE *fp) {
	int vec_len = vec.rows();
	for (int i = 0; i < vec_len; i++) {
		fprintf(fp, "%.15f\n", vec(i));
	}
}
