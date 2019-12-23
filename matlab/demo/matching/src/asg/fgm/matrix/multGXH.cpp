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
	m = (int) Ind0[2 * len];
        n = (int) Ind0[2 * len + 1];
	
	Ind = new int[2 * len];
	for (int p = 0; p < len; ++p) {
	    Ind[p * 2] = (int) Ind0[p * 2] - 1;
	    Ind[p * 2 + 1] = (int) Ind0[p * 2 + 1] - 1;
	}
    }

    return Ind;
}

/* 
 * function Y = multGXH(IndG, X, IndH)
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
    int mG, nG, mH, nH;
    
    // IndG
    double *IndG0 = mxGetPr(prhs[0]);
    int lenG = mxGetN(prhs[0]) - 1;
    int* IndG = indDou2Int(IndG0, lenG, m, mG, nG);

    // IndH
    double *IndH0 = mxGetPr(prhs[2]);
    int lenH = mxGetN(prhs[2]) - 1;
    int* IndH = indDou2Int(IndH0, lenH, n, mH, nH);

    // printf("mG %d nG %d mH %d nH %d m %d n %d\n", mG, nG, mH, nH, m, n);

    // check the dimension
    if (m != nG || n != mH) {
	mexErrMsgTxt("incorrect dimension");
    }

    // Y
    plhs[0] = mxCreateDoubleMatrix(mG, nH, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    int iG, jG, iH, jH, idxY, idxX;
    for (int pG = 0; pG < lenG; ++pG) {
	iG = IndG[2 * pG];
	jG = IndG[2 * pG + 1];
	
	for (int pH = 0; pH < lenH; ++pH) {
	    iH = IndH[2 * pH];
	    jH = IndH[2 * pH + 1];

	    idxY = jH * mG + iG;
	    idxX = iH * nG + jG;
	    Y[idxY] += X[idxX];
	}
    }

    // release memory
    delete[] IndG;
    delete[] IndH;
}
