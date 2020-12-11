#ifndef CPP_READ_MATRIX_UTILS_H
#define CPP_READ_MATRIX_UTILS_H

#include "LPboxADMMsolver.h"

/* The format of the dense matrix should be
 * 1,2,3
 * 6,7,8
 * for a 2x3 matrix
 */
void readDenseMat(FILE *fp, DenseMatrix &mat, int rowNum, int colNum);

/* Reading the dense matrix with the extra step of transforming the dense matrix
 * to sparse matrix
 */
void readDenseMatAsSparse(FILE *fp, SparseMatrix &smat, int rowNum, int colNum);

/* The format of the dense vector should be
 * 1
 * 2
 * 3
 * for a 3x1 vector
 */
void readDenseVec(FILE *fp, DenseVector &vec, int vecLen);

/* The format of the sparse matrix should be 
 * 1,1,1.5
 * 1,2,1.3
 * where the first index is the row index and the second index 
 * is the column index
 */
void readSparseMat(FILE *fp, SparseMatrix &mat);

/* Read the sparse matrix with diagonal elements initialized */
void readSparseMat_diag(FILE *fp, SparseMatrix &mat);

/* Generate the constraint matrices in the clustering task. The formulation of A, b, C, d can be 
 * found below equation (28)
 */
void generateMatrixFromLaplacian(const DenseMatrix &laplacian, SparseMatrix &A, DenseVector &b,
        SparseMatrix &C, DenseVector &d, int numNodes, int numClasses);

/* Writing a vector to the file */
void writeVectorToFile(const DenseVector &vec, FILE *fp);


#endif //CPP_READ_MATRIX_UTILS_H
