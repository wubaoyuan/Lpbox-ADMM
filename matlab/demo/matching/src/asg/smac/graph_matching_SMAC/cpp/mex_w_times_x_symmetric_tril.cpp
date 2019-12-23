/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23

mex_w_times_x_symmetric.c = used by ncuthard2.m in eigensolver.
 Examples:
    mex_w_times_x_symmetric(x,tril(A)) = A*x;
    A is sparse and symmetric, but x is full
Quicker with tril than with triu
Timothee Cour, Oct 12, 2003.
 
% test sequence:
    n=50;
    x=rand(n,1);
    A=sprand(n,n,0.01);
    A = A + A';
    y1 = A*x;
    y2 = mex_w_times_x_symmetric(x,tril(A));
    z = y1 - y2;
    norm(z(:))
 *=================================================================*/

# include "math.h"
# include "mex.h"
# include "a_times_b_cmplx.cpp"
/*# include "a_times_b.c"*/


void mexFunction(
int nargout,
mxArray *out[],
int nargin,
const mxArray *in[]
)
{
    int n;
    mwIndex *ir, *jc;
    double *x, *y, *pr;
        
    x = mxGetPr(in[0]);
    pr = mxGetPr(in[1]);
    ir = mxGetIr(in[1]);
    jc = mxGetJc(in[1]);
    
    n = mxGetM(in[1]);
    
    out[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    y = mxGetPr(out[0]);
    
    CSRsymm_VecMult_CAB_double_tril(n,pr,ir,jc,x,y);
}
