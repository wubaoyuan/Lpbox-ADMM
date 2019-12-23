#include "mex.h"

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define argmin(a, b)  (((a) < (b)) ? 0 : 1)

int* indDou2Int(double* Ind0, int& len, int k, int& m, int& n) {
    int* Ind;
    
    // empty
    if (len < 0) {
	m = k;
	n = k;
	len = k;
	
	Ind = new int[2 * len];
	for (int p = 0; p < len; ++p) {
	    Ind[p * 2] = p;
	    Ind[p * 2 + 1] = p;
	}
	// printf("len %d %d %d %d %d %d %d\n", len, Ind[0], Ind[1], Ind[2], Ind[3], Ind[38], Ind[39]);

    // non-empty
    } else {
	m = (int) Ind0[len * 2];
        n = (int) Ind0[len * 2 + 1];
	
	Ind = new int[len * 2];
	for (int p = 0; p < len; ++p) {
	    Ind[p * 2] = (int) Ind0[p * 2] - 1;
	    Ind[p * 2 + 1] = (int) Ind0[p * 2 + 1] - 1;
	}
    }

    return Ind;
}

/* 
 * function val = multGXHSQTr(indG, X, indH, IndS, Q)
 *
 *   X  -  m x n
 *   Y  -  mS x nS
 *   S  -  mS x nS
 *   G  -  mG x nG
 *   H  -  mH x nH
 *
 * History
 *   create  -  Feng Zhou, 03-20-2009
 *   modify  -  Feng Zhou, 09-03-2010
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // X
    double *X = mxGetPr(prhs[1]);
    int m = mxGetM(prhs[1]);
    int n = mxGetN(prhs[1]);
    
    // indG
    double *indG = mxGetPr(prhs[0]);
    int mG = mxGetM(prhs[0]);

    // indH
    double *indH = mxGetPr(prhs[2]);
    int nH = mxGetN(prhs[2]);

    // IndS
    double *IndS0 = mxGetPr(prhs[3]);
    int lenS = mxGetN(prhs[3]) - 1;
    int mS, nS;
    int* IndS = indDou2Int(IndS0, lenS, n, mS, nS);

    // Q
    double *Q = mxGetPr(prhs[4]);

    // printf("mG %d nG %d mH %d nH %d m %d n %d\n", mG, nG, mH, nH, m, n);

    // check the dimension
    if (mS != mG || nS != nH) {
	mexErrMsgTxt("incorrect dimension");
    }

    // val
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *val = mxGetPr(plhs[0]);
    *val = 0;

    int pS, iS, jS, i, j, idxX, idxY;
    for (int pS = 0; pS < lenS; ++pS) {
	iS = IndS[(pS << 1)];
	jS = IndS[(pS << 1) + 1];

	i = (int) indG[iS];
	i--;
	j = (int) indH[jS];
	j--;

	if (i < 0 || j < 0) {
	    continue;
	}

	idxY = jS * mS + iS;
	idxX = j * n + i;
	*val += X[idxX] * Q[idxY];
    }

    // release memory
    delete[] IndS;
}
