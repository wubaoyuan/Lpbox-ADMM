// function [X,scores] = mex_normalize_bistochastic(X,tol,maxIters);
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

#undef NDEBUG
#include <math.h>
#include <mex.h>
#include "mex_util.cpp"
#include "Matrix.cpp"

class Class1{
public:
    Matrix<double> X1;
    Matrix<double> X2;
    Matrix<double> Drowsum;
    Matrix<double> Dcolsum;
    Matrix<double> scores;
    double tol;
    int maxIters;
    int p, q;
    double eps;
    
    void mexGate(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
        X1.readmxArray(in[0]);
        p=X1.p;
        q=X1.q;
        eps=mxGetEps();
        tol=mxGetScalar(in[1]);
        maxIters=(int)mxGetScalar(in[2]);
        normalize_bistochastic();
        out[0]=X1.tomxArray();
        out[1]=scores.tomxArray();
	//         out[1]=mxCreateDoubleScalar((double)nbIters);
    }
    
    void normalize_bistochastic() {
        X2.resize(p,q);
        Drowsum.resize(p,1);
        Dcolsum.resize(q,1);
        int i;
        vector<double>scores0;
        for(i=0;i<maxIters;i++){
            copy_X1_to_X2(X1,X2);
            normalize_rows();
            normalize_columns();
            double score=compute_difference(X1,X2);
            scores0.push_back(score);
            if(score<tol){
                break;
            }
        }
        scores.read_vector(scores0);        
    }
    
    void normalize_rows(){
        compute_rowsum(X1);
        invert_D(Drowsum);
        int k=0;
        for(int j=0;j<q;j++){
            for(int i=0;i<p;i++,k++){
                X1[k]*=Drowsum[i];
            }
        }                
    }
    
    void normalize_columns(){
        compute_colsum(X1);
        invert_D(Dcolsum);
        int k=0;
        for(int j=0;j<q;j++){
            double temp=Dcolsum[j];
            for(int i=0;i<p;i++,k++){
                X1[k]*=temp;
            }
        }                                
    }
    double compute_difference(Matrix<double>&Y1,Matrix<double>&Y2){
        double score=0;
        int n=p*q;
        double temp=0;
        for(int i=0;i<n;i++){
            temp=Y1[i]-Y2[i];
            score+=temp*temp;
        }
        score=sqrt(score);
        return score;
    }
    void compute_rowsum(Matrix<double>&Y){
        Drowsum.set(0);
        int k=0;
        for(int j=0;j<q;j++){
            for(int i=0;i<p;i++,k++){
                Drowsum[i]+=Y[k];
            }
        }                
    }
    void compute_colsum(Matrix<double>&Y){
        int k=0;
        for(int j=0;j<q;j++){
            double temp=0;
            for(int i=0;i<p;i++,k++){
                temp+=Y[k];
            }
            Dcolsum[j]=temp;            
        }                                
    }  
    void invert_D(Matrix<double>&D){
        int n=D.n;
        for(int i=0;i<n;i++){
	    //             D[i]=1/(D[i]+eps);
            if(D[i])
                D[i]=1/D[i];
        }
    }
    void copy_X1_to_X2(Matrix<double>&Y1,Matrix<double>&Y2){
        int n=Y1.n;
        for(int i=0;i<n;i++)
            Y2[i]=Y1[i];
    }
};

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    Class1 class1;
    class1.mexGate(nargout,out,nargin,in);
}
