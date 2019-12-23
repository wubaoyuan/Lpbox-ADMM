/*
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

D = mex_computeRowSum(W);
D = sum(W,2);
handles both lower triangular or symmetric 
if full, assumes symmetric (not full)
 */

#include <mex.h>
#include <math.h>

#include "mex_util.cpp"

void computeColumnSum(double *rowSum,mwIndex* ir,mwIndex* jc,double*pr,int n){
    int end,i,j;
    for (j=0; j!=n; j++) {
        end = jc[j+1];
        for (i=jc[j]; i!=end; i++)
            *rowSum += *pr++;
        rowSum++;
    }
    rowSum -=n;
}
void computeRowSum(double *columnSum,mwIndex* ir,mwIndex* jc,double*pr,int n){
    int end,i,j;
    for (j=0; j!=n; j++) {
        end = jc[j+1];
        for (i=jc[j]; i!=end; i++)
            columnSum[ir[i]] += *pr++;
    }
}

void computeRowSum_tril(double *rowSum,mwIndex* ir,mwIndex* jc,double*pr,int n){
    int end,i,j;
    for (j=0; j!=n; j++) {
        i = jc[j];
        end = jc[j+1];
        for (; i<end && ir[i]<j; i++){
            pr++;
        }
		if(i<end && ir[i]==j){
            rowSum[j] += *pr++;
			i++;
		}
        for (; i!=end; i++) {
            rowSum[j] += *pr;
            rowSum[ir[i]] += *pr++;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n;
    mwIndex* ir;
    mwIndex *jc;
    double *pr,*D;
    
    if (mxIsSparse(prhs[0])) {
        n = mxGetM(prhs[0]);
        pr = mxGetPr(prhs[0]);
        ir = mxGetIr(prhs[0]);
        jc = mxGetJc(prhs[0]);
        
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        D = mxGetPr(plhs[0]);
        //computeRowSum(D,ir,jc,pr,n);
        //computeRowSum(D,ir,jc,pr,n);
        computeRowSum_tril(D, ir, jc, pr, n);

    } else {
            mxArray *args[2];
            args[0] = (mxArray*) prhs[0];
            args[1] = mxCreateDoubleScalar(2.0);
            mexCallMATLAB(1, plhs, 2, args, "sum");
        }
    }
    
