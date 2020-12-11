#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "solver.h"

const int n = 2;
int main() {
    std::ifstream Asource;
    std::ifstream bsource;
    std::ifstream xsource;

    Asource.open("(5041, 5041)", std::ios_base::in);
    bsource.open("b", std::ios_base::in);
    xsource.open("x0", std::ios_base::in);

//    for (int i = 0; i < 5041; i++) {
//        x_init[i] = ((double) rand() / RAND_MAX);
//        std::cout << x_init[i] << std::endl;
//    }

    if (Asource.is_open()) {
        std::cout << "The file is open" << std::endl;
    } else {
        std::cout << "The file is not open" << std::endl;
    }

    if (bsource.is_open()) {
        std::cout << "The file is open" << std::endl;
    } else {
        std::cout << "The file is not open" << std::endl;
    }
    int i = 0;

//    double *A = new double[5041 * 5041];
//    double *b = new double[5041];
//    double *x_init = new double[5041];
//
//    for (std::string line; std::getline(Asource, line);) {
//        if (line.empty()) {
//            break;
//        }
//        double d = std::stod(line);
//        A[i] = d / 2;
//        i++;
//    }
//
//    i = 0;
//    for (std::string line; std::getline(bsource, line);) {
//        if (line.empty()) {
//            break;
//        }
//        double d = std::stod(line);
//        b[i] = d;
//        i++;
//    }

    auto *A = new double[5041 * 5041];
    auto *b = new double[5041];
    auto *x_init = new double[5041];

    for (std::string line; std::getline(Asource, line);) {
        if (line.empty()) {
            break;
        }
        double d = std::stod(line);
        A[i] = d / 2;
        i++;
    }

    i = 0;
    for (std::string line; std::getline(bsource, line);) {
        if (line.empty()) {
            break;
        }
        double d = std::stod(line);
        b[i] = d;
        i++;
    }

    i = 0;
    for (std::string line; std::getline(xsource, line);) {
        if (line.empty()) {
            break;
        }
        double x = std::stod(line);
        x_init[i] = x;
        i++;
    }

//    int iteration_num = 5;
//    if (argc >= 2) {/
//        iteration_num = atoi(argv[1]);
//    }
    solver s;
    s.ADMM_bqp_unconstrained_init();
    s.set_x0(x_init);
//    double A[n * n] = {-0.5, 0, 0, -0.5};
//    double b[n] = {-1, -1};
    Solution sol;
//    s.set_max_iters(iteration_num);
//    double x0[2] = {1, 0};
    s.ADMM_bqp_unconstrained(0, 5041, A, b, sol, true);
    std::cout << sol.x_sol << std::endl;
    return 0;
}