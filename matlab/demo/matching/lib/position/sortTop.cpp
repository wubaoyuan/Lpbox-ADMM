#include "mex.h"
#define SWAP(a, b) {double temp = a; a = b; b = temp;}

static void minheap_pushdown(double *h, double *idx, int hs, int i);
static void minheap_raise(double *h, double *idx, int i);
void topK(double *top_k, double *top_idx, int k, const double *l, int n);

/*
 * function [bs, idx] = sortTop(as, m)
 * 
 * Get the position of the top m elements.
 *
 * Input
 *   as      -  array, 1 x n
 *   m       -  #elements
 *
 * Output
 *   bs      -  value, 1 x m
 *   idx     -  position, 1 x m
 *
 * History
 *   create  -  Feng Zhou, 02-28-2010
 *   modify  -  Feng Zhou, 02-28-2010
 */
void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {

    // as
    double* as = mxGetPr(prhs[0]);
    
    // m
    int m = int(*mxGetPr(prhs[1]));

    // bs
    plhs[0] = mxCreateDoubleMatrix(1, m, mxREAL);
    double *bs = mxGetPr(plhs[0]);

    // idx
    plhs[1] = mxCreateDoubleMatrix(1, m, mxREAL);
    double *idx = mxGetPr(plhs[1]);
    
    // n
    int n = mxGetN(prhs[0]);

    // algorithm
    topK(bs, idx, m, as, n);

    // reverse
    double tmp;
    int j = m - 1;
    for (int i = 0; i < j; ++i) {
	tmp = bs[i];
	bs[i] = bs[j];
	bs[j] = tmp;

	tmp = idx[i];
	idx[i] = idx[j];
	idx[j] = tmp;

	--j;
    }
}

static void minheap_pushdown(double *h, double *idx, int hs, int i) {
    int j = 0;

    if (2 * i + 2 < hs) {
	j = (h[2 * i + 1] < h[2 * i + 2]) ? 2 * i + 1 : 2 * i + 2;
    } else if (2 * i + 1 < hs) {
	j = 2 * i + 1;
    }

    if (j != 0 && h[j] < h[i]) {
	SWAP(h[i], h[j]);
	SWAP(idx[i], idx[j]);
	minheap_pushdown(h, idx, hs, j);
    }
}
 
static void minheap_raise(double *h, double *idx, int i) {

    if (i == 0) {
	return;
    }

    if (h[i] < h[(i - 1) / 2]) {
	SWAP(h[i], h[(i - 1) / 2]);
	SWAP(idx[i], idx[(i - 1) / 2]);
	minheap_raise(h, idx, (i - 1) / 2);
    }
}

void topK(double *top_k, double *top_idx, int k, const double *l, int n) {
    int i;

    for (i = 0; i < k; i++) {
	top_k[i] = l[i];
	top_idx[i] = i + 1;
	minheap_raise(top_k, top_idx, i);
    }

    for (i = k; i < n; i++) {
	if (l[i] > top_k[0]) {
	    top_k[0] = l[i];
	    top_idx[0] = i + 1;
	    minheap_pushdown(top_k, top_idx, k, 0);
	}
    }
}
