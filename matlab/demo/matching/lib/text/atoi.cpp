#include "mex.h"
#include "string.h"

#define MAX_LEN 100

/*
 * function i = atoi(a)
 *
 * History
 *   create  -  Feng Zhou, 08-15-2009
 *   modify  -  Feng Zhou, 08-15-2009
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // a
    char a[MAX_LEN];
    mxGetString(prhs[0], a, MAX_LEN);

    // main parts
    int i = atoi(a);

    // i
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = i;
}
