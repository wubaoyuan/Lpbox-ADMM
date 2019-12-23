/*
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

function istril = mex_istril(W);
istril = W==tril(W);
W : sparse matrix
 */


#include <mex.h>
#include <math.h>
#include "mex_util.cpp"

int aux(mwIndex* ir,mwIndex* jc,int n){
    int j,k,last;
    for (j=0; j!=n; j++) {
        last = jc[j+1];
        for (k=jc[j]; k!=last; k++) 
            if (ir[k] < j) 
                return 0;        
    }
    return 1;    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    assert(mxIsSparse(prhs[0]));
    int n;
    mwIndex *ir, *jc;
    
    n = mxGetN(prhs[0]);
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    
    plhs[0] = mxCreateDoubleScalar((double) aux(ir,jc,n));
}
