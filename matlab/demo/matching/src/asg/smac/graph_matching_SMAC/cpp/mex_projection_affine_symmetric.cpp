/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

 * mex_projection_affine_symmetric.cpp = used by computeEigenvectorsMRF in eigensolver.
 * mex_projection_affine_symmetric(X,W,constraint) = K W K X
 * L1 = constraint.L1;
 * L2 = constraint.L2;
 * B = constraint.B;
 * Ainv = constraint.Ainv;
 * Constraint : Cx=b
 * Ceq=[C(1:k-1,:)-(1/b(k))*b(1:k-1)*C(k,:)];
 * K = I-Ceq'*inv(Ceq*Ceq')*Ceq;
 *   = I-B'*Ainv*B + L1*L2';
 * Sherman-Morrison Formula: (http://mathworld.wolfram.com/Sherman-MorrisonFormula.html)
 *
 * W, B, Ainv are sparse, but L1, L2 are full
 * W is symmetric
 *
 * computation :
 
 *=================================================================*/

# include "math.h"
# include "mex.h"
#include <string.h> /* needed for memcpy() */

# include "mex_util.cpp"
# include "a_times_b_cmplx.cpp"
# include "mex_math.cpp"
/*# include "a_times_b.c"*/

// map<string,long>timing;
// long timei;
    
void compute_Kx(double *x, double *yn1, double *yn2, double *ym1, double *ym2, const mxArray *B,const mxArray *Ainv, double *L1,double *L2,int m, int n)
{
    //computes x:=x-B'*Ainv*B*x -L1*L2'*x
    
    //yn1=-(L2(:,1)'*x)*L1(:,1);
    scalar_times_vec(-dot(L2, x, n), L1, yn1, n);
    //yn2=-(L2(:,2)'*x)*L1(:,2);    
    scalar_times_vec(-dot(L2 + n, x, n), L1 + n, yn2, n);
    //yn1=yn1+yn2;
    vec_plus_vec(yn1, yn2, yn1, n);
    
    CSC_VecMult_CAB_double(m, n, mxGetPr(B), mxGetIr(B), mxGetJc(B), x, ym1);
    CSRsymm_VecMult_CAB_double(m, m, mxGetPr(Ainv), mxGetIr(Ainv), mxGetJc(Ainv), ym1, ym2);
    CSR_VecMult_CaABC_double(m, n, -1, mxGetPr(B), mxGetIr(B), mxGetJc(B), ym2, yn1);
    //yn1=x+yn1;
    vec_plus_vec(x,yn1,x,n);
    
}
void compute_KWKx(double *x, double *y, double *yn1, double *yn2, double *ym1, double *ym2, const mxArray *W, const mxArray *B, const mxArray *Ainv, double *L1,double *L2, int m, int n)
{
//     timei=clock();
    compute_Kx(x, yn1, yn2, ym1, ym2, B, Ainv, L1, L2, m, n);
//     timing["1"]+=clock()-timei;
//     timei=clock();

    CSRsymm_VecMult_CAB_double(n, n, mxGetPr(W), mxGetIr(W), mxGetJc(W), x, y);
//     timing["2"]+=clock()-timei;
//     timei=clock();

    compute_Kx(y, yn1, yn2, ym1, ym2, B, Ainv, L1, L2, m, n);
//     timing["1"]+=clock()-timei;
}

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[])
{
    //double *x = mxGetPr(mxDuplicateArray(in[0]));
    const mxArray *W = in[1];

    //const mxArray *Ceq = in[2];
    const mxArray *B = mxGetField(in[2], 0, "B");

    //const mxArray *Z = in[3];
    const mxArray *Ainv = mxGetField(in[2], 0, "Ainv");

    //double *b = mxGetPr(in[4]);
    double *L1 = mxGetPr(mxGetField(in[2], 0, "L1"));

    //double *a = mxGetPr(in[5]);
    double *L2 = mxGetPr(mxGetField(in[2], 0, "L2"));
    
    int m = mxGetM(B);
    int n = mxGetN(B);
    
    //double *x = (double*) memcpy(mxCalloc(n, sizeof(double)), mxGetPr(in[0]), sizeof(double)*n);
    double *x = copyArray(in[0]);
    
    out[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *y = mxGetPr(out[0]);
    double *yn1 = (double*) mxCalloc(n, sizeof(double));
    double *yn2 = (double*) mxCalloc(n, sizeof(double));
    double *ym1 = (double*) mxCalloc(m, sizeof(double));
    double *ym2 = (double*) mxCalloc(m, sizeof(double));
    
    compute_KWKx(x, y, yn1, yn2, ym1, ym2, W, B, Ainv, L1, L2, m, n);

    mxFree(yn1);
    mxFree(yn2);
    mxFree(ym1);
    mxFree(ym2);
    mxFree(x);
}
