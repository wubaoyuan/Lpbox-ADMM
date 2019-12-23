#include "mex.h"
#include "string.h"

#define MAX_LEN 100000
#define MAX_NUM_TOKEN 10000

/*
 * function tokens = tokenise(s, delim)
 *
 * History
 *   create  -  Feng Zhou, 08-15-2009
 *   modify  -  Feng Zhou, 08-15-2009
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // s
    char s[MAX_LEN];
    mxGetString(prhs[0], s, MAX_LEN);

    // delim
    char delim[MAX_LEN];
    mxGetString(prhs[1], delim, MAX_LEN);

    // tokens storage
    char* tokens[MAX_NUM_TOKEN];
    int m = 0;

    // main part
    char *token;

    token = strtok(s, delim);
    while (token != NULL) {
        tokens[m++] = token;
        token = strtok(NULL, delim);
    }

    // tokens
    plhs[0] = mxCreateCellMatrix(1, m);
    for (int i = 0; i < m; ++i) {
        mxSetCell(plhs[0], i, mxCreateString(tokens[i]));
    }
}
