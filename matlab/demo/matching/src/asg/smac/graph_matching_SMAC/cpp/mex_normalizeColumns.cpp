/*
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

mex_normalizeColumns(X);
in-place normalization of columns of X
[n,k] = size(X);
for  i=1:k
    X(:,i) = (X(:,i) / norm(X(:,i)));
end
 */

#include <mex.h>
#include <math.h>
#include "mex_util.cpp"
#include "Matrix.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n,k;
    int i,j,c;
    double sum;
    double *X;
    double *Y;
    
    n = mxGetM(prhs[0]);
    k = mxGetN(prhs[0]);

//     X = mxGetPr(prhs[0]);
    X = Matrix<double>::getData(prhs[0]);
    
    c = 0;
    
    if (nlhs == 0)
        
        for (j=0; j<k; j++) {
            sum = 0;
            for (i=0; i!=n; i++)
                sum += X[c] * X[c++];
            c -= n;
            sum = 1.0/sqrt(sum);
            //mexPrintf("sum = %6.5f\n",sum);
            for (i=0; i!=n; i++)
                X[c++] *= sum;
        }
    else {
        plhs[0] = mxCreateDoubleMatrix(n,k,mxREAL);
        Y = mxGetPr(plhs[0]);

        for (j=0; j<k; j++) {
            sum = 0;
            for (i=0; i!=n; i++)
                sum += X[c] * X[c++];
            c -= n;
            sum = 1.0/sqrt(sum);
            //mexPrintf("sum = %6.5f\n",sum);
            for (i=0; i!=n; i++)
                Y[c] = X[c++] * sum;
        }
    }
}

