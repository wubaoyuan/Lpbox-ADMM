#define OPENCV_DISABLE_EIGEN_TENSOR_SUPPORT

#include <unistd.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "image_segmentation_utils.h"
#include <getopt.h>

using namespace cv;


int main(int argc, char *argv[]) {
	int o;
	char* imagePath;
	int numNodes = 5000;
	char* outputPath;
    std::string logFile;
    bool log_file_is_set = false;
	int option_idx;

	struct option long_options[] = {
        {"input", required_argument, NULL, 'i'},
        {"pixel", required_argument, NULL, 'n'},
        {"output", required_argument, NULL, 'o'},
        {"log", required_argument, NULL, 'l'},
    };

	if (argc == 1) {
		printf("Usage: image_segmentation <command>\n\nCommands:\n  -i --input <path>\tInput file path\n"
				"  -n --pixel <num_pixels>\tNumber of pixels of the output picture (5000 by default)\n"
				"  -o --output <path>\tOutput file path (prefer .bmp extension to avoid opencv bug)\n"
				"  -l --log <path>\tLog file path \n");
		return 0;
	}
    const char *optstring = "i:n:o:l:h";
    while ((o = getopt_long(argc, argv, optstring, long_options, &option_idx)) != -1) {
        switch (o) {
            case 'i':
				imagePath = optarg;
                break;
            case 'n':
                numNodes = atoi(optarg);
				break;
            case 'o':
				outputPath = optarg;
				break;
			case 'l':
				log_file_is_set = true;
				logFile = optarg;
				break;
			default:
				printf("Usage: image_segmentation <command>\n\nCommands:\n  -i --input <path>\tInput file path\n"
						"  -n --pixel <num_pixels>\tNumber of pixels of the output picture (5000 by default)\n"
						"  -o --output <path>\tOutput file path (prefer .bmp extension to avoid opencv bug)\n"
						"  -l --log <path>\tLog file path \n");
				return 0;
		}
	}
	Mat image = imread(imagePath, 0);
	double lambda = 3.0;

	/* Scaling the picture */
	double scale = std::sqrt(numNodes / (double) (image.rows * image.cols));
	Mat scaled_image;
	Size scaled_size = Size(std::round(scale * image.cols), std::round(scale * image.rows));
	resize(image, scaled_image, Size(), scale, scale);

	int numPixels = scaled_image.rows * scaled_image.cols;
	DenseMatrix I;
	cv2eigen(scaled_image, I);

	/* Rescale the value of each pixel between [0, 1] */
	I.array() = I.array() / 256.0;
	DenseMatrix unary_cost;
	SparseMatrix binary_cost;

	double sigma = 0.1;
	double b = 0.6;
	double f1 = 0.2;
	double f2 = 0.2;

	get_unary_cost(I, sigma, b, f1, f2, unary_cost);
	get_binary_cost(I, binary_cost); /* Actually calculating round(lambda * W) in this so we can get the 
									  * binary cost directly */

	unary_cost.array() = unary_cost.array().round();

	SparseMatrix _A;
	DenseVector _b;
	get_A_b_from_cost(unary_cost, binary_cost, _A, _b);

	printf("Finished generating matrices, starting algorithm\n");
	Solution sol;
	DenseVector x0 = DenseVector::Zero(_b.rows());
	LPboxADMMsolver solver;
	solver.ADMM_bqp_unconstrained_init();

	if (log_file_is_set) {
		solver.set_log_file(logFile);
		printf("Writing log output to %s\n", logFile.c_str());
	}
	solver.ADMM_bqp_unconstrained_legacy(_A.cols(), _A, _b, x0, sol);
	printf("Algorithm finishes, saving picture\n");

	/* Reshape the solution to fit the size of the image */
	DenseMatrix reshaped_mat = Eigen::Map<DenseMatrix>(sol.x_sol->data(), I.rows(), I.cols());

	/* Setting the location with value 1 to be white */
	DenseIntMatrix r_mat = (reshaped_mat.array() >= 0.5).matrix().cast<int>() * 255;
	Mat res_mat(r_mat.rows(), r_mat.cols(), CV_8UC1);
	eigen2cv(r_mat, res_mat);
	imwrite(outputPath, res_mat);

	delete sol.best_sol;
	delete sol.x_sol;
	delete sol.y1;
	delete sol.y2;
	return 1;
}
