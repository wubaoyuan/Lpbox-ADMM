#include "mex.h"
#include "stdlib.h"

#define MAX_LEN 100

/*
 * function f = atof(a)
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
    double f = atof(a);

    // f
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = f;
}
