//
// Created by squintingsmile on 2020/10/5.
//

#include <iostream>
#include <unistd.h>
#include <string>
#include "solver.h"

struct Trip {
    int row;
    int col;
    float_t val;
};

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Please enter the read and send pipe file descriptor" << std::endl;
    }

    int readPipe = atoi(argv[1]);
    int sendPipe = atoi(argv[2]);
    int bLen = 0;
    auto dataLen = read(readPipe, &bLen, sizeof(int));
    if (dataLen < 0) {
        perror("read()");
        exit(-1);
    }
    std::cout << "bLen: " << bLen << std::endl;

    float_t *x0 = new float_t[bLen];
    dataLen = read(readPipe, x0, sizeof(float_t) * bLen);
    if (dataLen < 0) {
        perror("read()");
        exit(-1);
    }
    std::cout << "received x0" << std::endl;

    float_t *b = new float_t[bLen];
    dataLen = read(readPipe, b, sizeof(float_t) * bLen);
    if (dataLen < 0) {
        perror("read()");
        exit(-1);
    }
    DenseVector _b(bLen);
    memcpy(_b.data(), b, sizeof(float_t) * bLen);
    std::cout << "received b" << std::endl;

    int ALen = 0;
    dataLen = read(readPipe, &ALen, sizeof(int));
    if (dataLen < 0) {
        perror("read()");
        exit(-1);
    }
    std::cout << "received ALen" << std::endl;

    std::vector<Triplet> tripArr;
    tripArr.reserve(ALen);
    Trip tr;
    std::cout << "ALen: " << ALen << std::endl;
    for (int i = 0; i < ALen; i++) {
//        std::cout << "i: " << i << std::endl;
        read(readPipe, &tr, sizeof(tr));
//        std::cout << "i: " << i << std::endl;
//        std::cout << "row: " << tr.row << std::endl;
//        std::cout << "col: " << tr.col << std::endl;
//        std::cout << "val: " << tr.val << std::endl;
        tripArr.push_back(Triplet(tr.row, tr.col, tr.val));
    }
    std::cout << "size of A: " << tripArr.size() << std::endl;
    SparseMatrix A(bLen, bLen);
    A.setFromTriplets(tripArr.begin(), tripArr.end());

    solver solv;
    Solution sol;
    solv.ADMM_bqp_unconstrained_init();
    solv.set_x0(x0);
    solv.ADMM_bqp_unconstrained(0, bLen, A, _b, sol, true);

    int res = write(sendPipe, sol.x_sol, sizeof(float_t) * bLen);
    if (res < 0) {
        perror("write()");
        exit(-1);
    }

    return 1;
}