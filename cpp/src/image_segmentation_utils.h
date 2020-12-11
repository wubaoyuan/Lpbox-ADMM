#ifndef CPP_IMAGE_SEGMENTATION_UTIL_H
#define CPP_IMAGE_SEGMENTATION_UTIL_H

#include "LPboxADMMsolver.h"

void get_unary_cost(const DenseMatrix &nodes, _double_t sigma, _double_t b, _double_t f1, _double_t f2, DenseMatrix &unary_cost);

void get_binary_cost(const DenseMatrix &image, SparseMatrix &binary_cost);

void get_A_b_from_cost(const DenseMatrix &unary_cost, const SparseMatrix &binary_cost, SparseMatrix &A, DenseVector &b);

#endif //CPP_IMAGE_SEGMENTATION_UTIL_H
