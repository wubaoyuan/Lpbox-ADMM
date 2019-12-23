#include "mex.h"

/*
 * function varargout = oness(varargin)
 *
 * Input
 *   varargin   -  dimensions
 *
 * Output
 *   varargout  -  cell matrices
 *
 * History
 *   create     -  Feng Zhou, 02-13-2009
 *   modify     -  Feng Zhou, 07-30-2009
 */
void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {

    for (int i = 0; i < nlhs; ++i) {
        mexCallMATLAB(1, plhs + i, nrhs, const_cast<mxArray**>(prhs), "ones");
    }
}
