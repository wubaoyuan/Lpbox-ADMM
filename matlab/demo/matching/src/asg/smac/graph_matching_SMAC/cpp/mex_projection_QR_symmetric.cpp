/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

 * mex_projection_QR_symmetric.cpp = used by quickNcutHardBiais.m in eigensolver.
 * mex_projection_QR_symmetric(X,P,Ubar,R) == (eye(length(P))-Ubar*(R'*R)^{-1}*Ubar') * P * (eye(length(P))-Ubar*(R'*R)^{-1}*Ubar') * X ;
 * (R'*R)^{-1} is solved by solving a triangular system
 
 * P, Ubar, R are sparse, but X is full
 
 * R is upper triangular
 * P is symmetric
 
% test sequence:
 
 *=================================================================*/

# include "math.h"
# include "mex.h"
# include "a_times_b_cmplx.cpp"
/*# include "a_times_b.c"*/

//#include "mex_util.cpp"
//long timei,time0[7];

void mexFunction(
int nargout,
mxArray *out[],
int nargin,
const mxArray *in[]
)
{
	//int i;
	//for(i=0;i<6;i++)
	//	time0[i]=0;

	//timei=clock();

    mwIndex *ir[4], *jc[4];
    int m[4], n[4];
    double *y, *y1,*y2,*y3,*y4,*y5,*y6, *pr[4];
    double *y2bis, *y5bis;
    int k;
    
    for (k=0; k<4; k++) {
        m[k] = mxGetM(in[k]);
        n[k] = mxGetN(in[k]);
        pr[k] = mxGetPr(in[k]);
        if (k>0) {
            ir[k] = mxGetIr(in[k]);
            jc[k] = mxGetJc(in[k]);
        }
    }
    
    out[0] = mxCreateDoubleMatrix(m[1],1,mxREAL);
    y = mxGetPr(out[0]);
    
    
    y1 = (double*)mxCalloc(n[2], sizeof(double));
    y2 = (double*)mxCalloc(m[3], sizeof(double));
    y2bis = (double*)mxCalloc(m[3], sizeof(double));
    y3 = (double*)mxCalloc(m[1], sizeof(double));
    y4 = (double*)mxCalloc(m[1], sizeof(double));
    y5 = (double*)mxCalloc(n[2], sizeof(double));
    y5bis = (double*)mxCalloc(n[2], sizeof(double));
    y6 = (double*)mxCalloc(n[2], sizeof(double));
    
    CSR_VecMult_CAB_double(m[2],n[2],pr[2],ir[2],jc[2],pr[0],y1);
    CSR_VecTriangSlvLD_CAB_double(m[3],pr[3],ir[3],jc[3],y1,y2bis);
    CSC_VecTriangSlvUD_CAB_double(m[3],pr[3],ir[3],jc[3],y2bis,y2);
    for(k=0;k<m[1];k++) {
        y3[k]=pr[0][k];
    }
    CSC_VecMult_CaABC_double(m[2],n[2],-1,pr[2],ir[2],jc[2],y2,y3);
    
    //CSRsymm_VecMult_CAB_double(m[1],n[1],pr[1],ir[1],jc[1],y3,y4);
	CSRsymm_VecMult_CAB_double_tril(m[1],pr[1],ir[1],jc[1],y3,y4);
	
	//time0[0]+=clock()-timei;
	//timei=clock();
	//CSRsymm_VecMult_CAB_double_tril(m[1],pr[1],ir[1],jc[1],y3,y4);
	//time0[1]+=clock()-timei;
  	//timei=clock();
  

    CSR_VecMult_CAB_double(m[2],n[2],pr[2],ir[2],jc[2],y4,y5);
    
    CSR_VecTriangSlvLD_CAB_double(m[3],pr[3],ir[3],jc[3],y5,y5bis);
    CSC_VecTriangSlvUD_CAB_double(m[3],pr[3],ir[3],jc[3],y5bis,y6);
    for(k=0;k<m[1];k++) {
        y[k]=y4[k];
    }
    CSC_VecMult_CaABC_double(m[2],n[2],-1,pr[2],ir[2],jc[2],y6,y);
    
    
    mxFree(y1);
    mxFree(y2);
    mxFree(y2bis);
    mxFree(y3);
    mxFree(y4);
    mxFree(y5);
    mxFree(y5bis);
    mxFree(y6);

	//time0[0]+=clock()-timei;
	//dispTime(time0,2);
}
