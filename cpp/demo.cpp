//
// Created by squintingsmile on 2020/10/5.
//

#include <unistd.h>
#include <sys/stat.h>
#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Please enter the input image" << std::endl;
    }

    char *inputImagePath = argv[1];
    char *outputImagePath;

    int ret = mkdir("out", 0777);
    if (ret == -1) {
        if (errno != EEXIST) {
            perror("mkdir()");
            return -1;
        }
    }

    if (argc < 3) {
        outputImagePath = new char[100];
        sprintf(outputImagePath, "/out/%s", inputImagePath);
    } else {
        outputImagePath = argv[2];
    }

    int pipefd1[2];
    int pipefd2[2];
    pipe(pipefd1);
    pipe(pipefd2);

    pid_t pid = fork();
    if (pid == -1) {
        perror("fork()");
        exit(EXIT_FAILURE);
    }

    if (pid == 0) {
        close(pipefd1[1]);
        close(pipefd2[0]);
        char command[100];
        sprintf(command, "python3 DataProcessing.py %d %d %s %s", pipefd1[0], pipefd2[1], inputImagePath, outputImagePath);
        system(command);
        close(pipefd1[0]);
        close(pipefd2[1]);
        exit(EXIT_SUCCESS);
    } else {
        close(pipefd1[0]);
        close(pipefd2[1]);
        char command[100];
        sprintf(command, "./cpp_algo %d %d", pipefd2[0], pipefd1[1]);
        system(command);
        close(pipefd2[0]);
        close(pipefd1[1]);
        exit(EXIT_SUCCESS);
    }

    return 1;
}
