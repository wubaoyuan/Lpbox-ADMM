#include "LPboxADMMsolver.h"
#include "read_matrix_utils.h"
#include <vector>
#include <unistd.h>
#include <iostream>
#include <getopt.h>

/* Each row is the corresponding vector. The index matrix's index starts from 0.
 * returns the nearest indices corresponds to each index and the corresponding distance
 */
void simpleKNN(const DenseMatrix &mat, int num_neighbors, DenseIntMatrix &idx_mat, DenseMatrix &dis_mat) {
	int numSamples = mat.rows();
	int dim = mat.cols();
	idx_mat = DenseIntMatrix::Zero(num_neighbors, numSamples);
	dis_mat = DenseMatrix::Zero(num_neighbors, numSamples);

	for (int i = 0; i < numSamples; i++) {
		std::vector<int> idx_vec(num_neighbors, 0);
		std::vector<double> dis_vec(num_neighbors, 100000.0);
		for (int j = 0; j < numSamples; j++ ) {
			if (i == j) {
				continue;
			}
			double distance = (mat.block(i, 0, 1, dim).matrix() - mat.block(j, 0, 1, dim).matrix()).squaredNorm();
			auto iter = std::max_element(dis_vec.begin(), dis_vec.end());
			if (distance < *iter) {
				int vec_idx = iter - dis_vec.begin();
				idx_vec[vec_idx] = j;
				*iter = distance;
			}
		}
		for (int idx = 0; idx < num_neighbors; idx++) {
			idx_mat(idx, i) = idx_vec[idx];
			dis_mat(idx, i) = dis_vec[idx];
		}
	}	
}


int main(int argc, char **argv) {
	
	int num_samples = 0;
	int num_classes = 0;
	int num_dim = 0;
	int num_neighbors = 0;
	char *inputFile;
	std::string outputFile;
	bool output_file_is_set = false;
	int permutation;

	int o;
	int option_idx;
	struct option long_options[] = {
		{"input", required_argument, NULL, 'i'},
		{"sample", required_argument, NULL, 's'},
		{"class", required_argument, NULL, 'c'},
		{"dim", required_argument, NULL, 'd'},
		{"neighbor", required_argument, NULL, 'n'},
		{"log", required_argument, NULL, 'l'},
		{"permute", required_argument, NULL, 'p'},
	};

    if (argc == 1) {
        printf("Usage: image_segmentation <command>\n\nCommands:\n  -i --input <path>\tInput file path\n"
                "  -s --sample <num_samples>\tNumber of samples\n"
                "  -c --class <num_classes>\tNumber of classes\n"
                "  -d --dim <num_dim>\tDimension of the data\n"
                "  -n --neighbor <neighbor>\tNumber of neighbors in KNN\n"
                "  -l --log <path>\tLog file path\n"
                "  -p --permute <permutation>\tpermutation of the data. Default is the index in the first column\n"
				"\t\tdata in the middle, and the label in the last column. 1 for label in the last \n"
				"\t\tcolumn and the data in the first few columns. 2 for label in the first column \n"
				"\t\tand data in the last few columns\n");
        return 0;
    }
    const char *optstring = "i:s:c:l:n:d:p:h";
    while ((o = getopt_long(argc, argv, optstring, long_options, &option_idx)) != -1) {
        switch (o) {
            case 'i':
                inputFile = optarg;
                break;
            case 's':
                num_samples = atoi(optarg);
                break;
            case 'c':
                num_classes = atoi(optarg);
                break;
            case 'n':
                num_neighbors = atoi(optarg);
                break;
            case 'd':
                num_dim = atoi(optarg);
                break;
            case 'l':
                outputFile = optarg;
				output_file_is_set = true;
                break;
            case 'p':
				permutation = atoi(optarg);
                break;
            default:
				printf("Usage: image_segmentation <command>\n\nCommands:\n  -i --input <path>\tInput file path\n"
						"  -s --sample <num_samples>\tNumber of samples\n"
						"  -c --class <num_classes>\tNumber of classes\n"
						"  -d --dim <num_dim>\tDimension of the data\n"
						"  -n --neighbor <neighbor>\tNumber of neighbors in KNN\n"
						"  -l --log <path>\tLog file path\n"
						"  -p --permute <permutation>\tpermutation of the data. Default is the index in the first column\n"
						"\t\tdata in the middle, and the label in the last column. 1 for label in the last \n"
						"\t\tcolumn and the data in the first few columns. 2 for label in the first column \n"
						"\t\tand data in the last few columns\n");
                return 0;
        }
    }

	FILE *input_mat_fp = fopen(inputFile, "r");
	
	DenseMatrix mat;

	/* For this task, the index is written in the first column and the ground truth labels are 
	 * written in the last column. If the data is permuted differently, please edit this line to
	 * fit your data
	 */
	DenseMatrix vecs;
	DenseIntMatrix indices;
	DenseIntMatrix true_labels;

	/* The first column is the index, the last column is the label. The middle is vector (for the glass dataset) */
	if (permutation == 0) {
		readDenseMat(input_mat_fp, mat, num_samples, num_dim + 2); 
		vecs = mat.block(0, 1, num_samples, num_dim);
		indices = (mat.block(0, 0, num_samples, 1).array().round()).cast<int>();
		true_labels = (mat.block(0, num_dim + 1, num_samples, 1).array().round()).cast<int>();
	}
	
	/* The first few columns are the vector, the last column is the label (for the iris dataset)*/
	if (permutation == 1) {
		readDenseMat(input_mat_fp, mat, num_samples, num_dim + 1);
		vecs = mat.block(0, 0, num_samples, num_dim);
		true_labels = (mat.block(0, num_dim, num_samples, 1).array().round()).cast<int>();
	}

	/* The first column is the label, followed by the vectors (for the wine dataset)*/
	if (permutation == 2) {
		readDenseMat(input_mat_fp, mat, num_samples, num_dim + 1);
		vecs = mat.block(0, 1, num_samples, num_dim);
		true_labels = (mat.block(0, 0, num_samples, 1).array().round()).cast<int>();
	}
	fclose(input_mat_fp);

	DenseIntMatrix idx_mat; /* Stores the index of the k nearest neighbor */
	DenseMatrix dis_mat; /* Stores the distance of the k nearest neighbor */


	/* Normalizing the feature values */
	for (int i = 0; i < num_dim; i++) {
		auto col = vecs.col(i);
		double max = col.maxCoeff();
		double min = col.minCoeff();
		double range = max - min;
		col.array() = 2 * (col.array() - min) / range - 0.5;

	}

	/* This KNN is a naive version and should be replaced by faster algorithm if the data is large */
	simpleKNN(vecs, num_neighbors, idx_mat, dis_mat);
	double max_dist = dis_mat.maxCoeff();
	DenseMatrix dis_mat_normalized = dis_mat.array() / max_dist;

	DenseMatrix dis_mat_squared = DenseMatrix::Zero(num_samples, num_samples);
	for (int i = 0; i < num_samples; i++) {
		for (int j = 0; j < num_neighbors; j++) {
			int idx = idx_mat(j, i);
			dis_mat_squared(idx, i) = dis_mat_normalized(j, i);
		}
	}

	dis_mat_squared = (dis_mat_squared + dis_mat_squared.transpose()) / 2;

	DenseMatrix W_matrix = DenseMatrix::Zero(num_samples, num_samples);
	for (int i = 0; i < dis_mat_squared.rows(); i++) {
		for (int j = 0; j < dis_mat_squared.cols(); j++) {
			if (dis_mat_squared(i, j) != 0.0) {
				W_matrix(i, j) = std::log2(dis_mat_squared(i, j));
			}
		}
	}
	 
	DenseMatrix M_matrix = - W_matrix;

	DenseMatrix D_matrix = DenseMatrix::Zero(num_samples, num_samples);
	D_matrix.diagonal().array() = (M_matrix.colwise().sum()).transpose();

	DenseMatrix laplacian = D_matrix - M_matrix;

	SparseMatrix A, C;
	DenseVector b, d;

	/* The actual constraint is induced from the laplacian matrix. Refer to the formula below (28) */
	generateMatrixFromLaplacian(laplacian, A, b, C, d, num_samples, num_classes);

	int b_len = b.rows();
	int d_len = d.rows();
	LPboxADMMsolver solv;
	Solution sol;
	DenseVector x0 = DenseVector(b_len);
	x0.setZero();
	solv.ADMM_bqp_linear_eq_init();
	solv.set_pcg_tol(1e-3);

	if (output_file_is_set) {
		solv.set_log_file(outputFile);
		printf("Writing log output to %s\n", outputFile.c_str());
	}

	solv.ADMM_bqp_linear_eq(b_len, A, b, x0, d_len, C, d, sol);

	/* Reshaping the x solution to a matrix with number of rows equal to num_samples and number of
	 * columns equal to num_classes
	 */
	DenseMatrix x_sol = Eigen::Map<DenseMatrix>((*sol.x_sol).data(), num_samples, num_classes);
	DenseMatrix res_labels = DenseMatrix::Zero(num_samples, num_classes);
	for (int i = 0; i < num_samples; i++) {
		int max_col;
		int max_row;
		x_sol.row(i).maxCoeff(&max_row, &max_col);
		res_labels(i, max_col) = 1.0; /* label starts from 0 */
	}
	
	double objective = (res_labels.transpose() * laplacian * res_labels).trace() - D_matrix.diagonal().sum();
	printf("objective is %lf\n", objective);

	delete sol.best_sol;
	delete sol.x_sol;
	delete sol.y1;
	delete sol.y2;
	
	return 1;
}
