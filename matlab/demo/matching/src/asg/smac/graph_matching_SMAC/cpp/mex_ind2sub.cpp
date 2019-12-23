//function [indi,indj]=mex_ind2sub([p,q],ind);
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.


# include "math.h"
# include "mex.h"

//# include "mex_math.cpp"
# include "mex_util.cpp"
# include "Matrix.cpp"

/*
typedef int Tkey;
void ind2sub(int*ind,int*indi,int*indj,int p,int q,int n){
    int i;
    int*pindi=indi;
    int*pindj=indj;
    int*pind=ind;
    int temp;
    for(i=0;i<n;i++){
        temp=*pind++-1;
 *pindi++=temp%p+1;
 *pindj++=temp/p+1;
    }
}
 */

/*
void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    Matrix<int>pq(in[0]);
    assert(pq.n==2);
    int p=pq[0];
    int q=pq[1];
    
    
    // 	double* pq=mxGetPr(in[0]);
    // 	int p=(int)pq[0];
    // 	int q=(int)pq[1];
    
    //     image=Matrix<T>::getData(in[0]);
    
    int*ind=readInt32(in[1]);
    int n=mxGetNumberOfElements(in[1]);
    out[0] = create_mxArray(n,1,mxINT32_CLASS);
    out[1] = create_mxArray(n,1,mxINT32_CLASS);
    int*indi=(int*)mxGetData(out[0]);
    int*indj=(int*)mxGetData(out[1]);
    ind2sub(ind,indi,indj,p,q,n);
}

*/


template<class T>class Class1{
    public:
        void ind2sub(T*ind,T*indi,T*indj,int p,int q,int n){
            //caution: output is in matlab conventions
            int i;
            T*pindi=indi;
            T*pindj=indj;
            T*pind=ind;
            T temp;
            for(i=0;i<n;i++){
                temp=*pind++-1;
                *pindi++=((int)temp)%p+1;
                *pindj++=((int)temp)/p+1;
            }
        }
        void mexGate(int nargout, mxArray *out[], int nargin, const mxArray *in[]){
            Matrix<int>pq(in[0]);
            assert(pq.n==2);
            int p=pq[0];
            int q=pq[1];
            
            T*ind=Matrix<T>::getData(in[1]);
            int n=mxGetNumberOfElements(in[1]);
            
            out[0]=Matrix<T>::createmxArray(n,1);
            out[1]=Matrix<T>::createmxArray(n,1);
            
            copyDims(in[1],out[0]);
            copyDims(in[1],out[1]);
            
            T*indi=Matrix<T>::getData(out[0]);
            T*indj=Matrix<T>::getData(out[1]);
            
            ind2sub(ind,indi,indj,p,q,n);
        }
        
};

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    mxClassID classID=mxGetClassID(in[1]);
    if(classID==mxDOUBLE_CLASS){
        Class1<double> class1; class1.mexGate(nargout,out,nargin,in);
    }
    else if(classID==mxSINGLE_CLASS){
        Class1<float> class1; class1.mexGate(nargout,out,nargin,in);
    }
    else if(classID==mxINT32_CLASS){
        Class1<int> class1; class1.mexGate(nargout,out,nargin,in);
    }
    else
        assert(0);
}


