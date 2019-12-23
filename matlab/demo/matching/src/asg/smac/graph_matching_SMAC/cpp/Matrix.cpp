// class Matrix
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.


#pragma once
#include <math.h>
#include <mex.h>
#include <algorithm>
#include <vector>
#include "mex_util.cpp"
using namespace std;

//namespace Matrix;
/*
 mxDOUBLE_CLASS

 TODO:
support better for bool (problem with set,add in *pdata++)
for ex, replace mxBOOL by char, etc
support for tomxArray() using _uninitialized... with caution
 */
template<class T> class Matrix{
public: 
	int p,q,n;
	vector<T>data;

	Matrix(){
		p=0;
		q=0;
		n=0;
		data.clear();
	}
	Matrix(int p,int q){
		this->p=p;
		this->q=q;
		n=p*q;
		data.resize(n);
	}
	Matrix(const mxArray *A){
		readmxArray(A);
	}
	~Matrix(){
	}
	T& operator[] (const int index){
		return data[index];
	}
	Matrix<T>& copy(){
		Matrix<T>matrix(p,q);
		//T*pdata=&data[0];
		//T*pdata2=&matrix.data[0];
		//int i=n;
		//while (--i >= 0) *pdata2++ = *pdata++;
		for(int i=0;i<n;i++)
			matrix.data[i] = data[i];
		return matrix;
	}
	void copy(Matrix<T>&matrix){
		resize(matrix.p,matrix.q);
		T*pdata=&data[0];
		T*pdata2=&matrix.data[0];
		int i=n;
		while (--i >= 0) *pdata++ = *pdata2++;
	}
	void read_vector(vector<T>&v){
		resize(v.size(),1);
		T*pdata=&data[0];
		T*pdata2=&v[0];
		int i=n;
		while (--i >= 0) *pdata++ = *pdata2++;
	}



	void resize(int p,int q){
		this->p=p;
		this->q=q;
		n=p*q;
		data.resize(n);
	}
	void push_back(T& val){
		n++;
		p++;//assumes column vector
		if(q==0){
			q=1;
		}
		data.push_back(val);
	}
	void erase(typename vector<T>::iterator iter){
		data.erase(iter);
		n--;
		p--;//assumes column vector
	}
	void add(T val){
		T*pdata=&data[0];
		int i=n;
		while (--i >= 0) *pdata++ += val;
		//while (--i >= 0) pdata[i] += val;
		/*
		if(val==1)
		while (--i >= 0) ++(*pdata++);
		else if(val==-1)
		while (--i >= 0) --(*pdata++);
		else
		while (--i >= 0) *pdata++ += val;
		*/
	}
	void set(T val){
		T*pdata=&data[0];
		int i=n;
		while (--i >= 0) *pdata++ = val;
	}
	void reset(){
		T*pdata=&data[0];
		int i=n;
		while (--i >= 0) *pdata++ = 0;
	}
	bool isSameSize(Matrix&matrix2){
		return matrix2.p==p && matrix2.q==q;
	}

	void readmxArray(const mxArray *A){
		//if A is pAxqAxkA, p=pA,q=qA*kA
		readmxArray(A,data);
		p=mxGetM2(A);
		q=mxGetN2(A);
		n=p*q;
	}
	static T* getData(const mxArray *A){
        if(mxIsEmpty(A)){
            T*vals=NULL;
            return vals;
        }
		assert(!mxIsSparse(A));
        mxClassID classID=mxGetClassID(A);
        assert(MatlabInterface::getClassID(typeid(T))==classID);
        T*vals=(T*)mxGetData(A);
        return vals;
    }
	static void readmxArray(const mxArray *A,vector<T>&data){
		assert(!mxIsSparse(A));
		int n=mxGetNumberOfElements(A);
		data.resize(n);
		int i=n;

		mxClassID classID=mxGetClassID(A);


		//vector<bool> special 
		if(classID == mxLOGICAL_CLASS){
			bool *pA=(bool*)mxGetData(A);
			for(i=0;i<n;i++)data[i] = *pA++;
			//BIG BUG !!
			//while (--i >= 0) data[i] = *pA++;
			return;
		}
		switch(classID){
		case mxDOUBLE_CLASS:{
			double *pA=mxGetPr(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
							}
		case mxSINGLE_CLASS:{
			float *pA=(float*)mxGetData(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
							}
		case mxINT32_CLASS:{
			int *pA=(int*)mxGetData(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
						   }
		case mxUINT32_CLASS:{
			unsigned int *pA=(unsigned int*)mxGetData(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
						   }
		case mxUINT16_CLASS:{
			unsigned short int *pA=(unsigned short int*)mxGetData(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
						   }
		case mxUINT8_CLASS:{
			unsigned char *pA=(unsigned char*)mxGetData(A);
			for(i=0;i<n;i++)
				data[i] = (T)*pA++;
			break;
						   }
		default: assert(0);
		}

	}

// 	static mxArray* tomxArray(vector<T>&v){
//         //TODO
//     }
	mxArray* tomxArray(){
		return tomxArray(MatlabInterface::getClassID(typeid(T)));
	}
	//template<>
	mxArray* tomxArray(mxClassID classID){

		//vector<bool> special 
		if(classID == mxLOGICAL_CLASS){
			mxArray *A=mxCreateLogicalMatrix(p,q);
			bool *pA=(bool*)mxGetData(A);

			//BIG BUG !!
			//int i=n;
			//while (--i >= 0) *pA++ = data[i];//BUG !! //does not work !!
			//while (--i >= 0) pA[i] = data[i];
			for(int i=0;i<n;i++)
				pA[i] = data[i];

			return A;
		}

		int ndim=2;
		const mwSize dims[]={(mwSize)p,(mwSize)q};
		mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);

		//T*pdata=&data[0];
		int i=n;
		switch(classID){
			case mxDOUBLE_CLASS:{
				double *pA=mxGetPr(A);
				for(i=0;i<n;i++) *pA++ = data[i];
				//while (--i >= 0) *pA++ = *pdata++;
				break;
								}
			case mxSINGLE_CLASS:{
				float *pA=(float*)mxGetData(A);
				for(i=0;i<n;i++) *pA++ = data[i];
				//while (--i >= 0) *pA++ = *pdata++;
				break;
								}
			case mxINT32_CLASS:{
				int *pA=(int*)mxGetData(A);
				for(i=0;i<n;i++) *pA++ = data[i];
				//while (--i >= 0) *pA++ = *pdata++;
				break;
							   }
			case mxUINT8_CLASS:{
				unsigned char *pA=(unsigned char*)mxGetData(A);
				for(i=0;i<n;i++) *pA++ = data[i];
				//while (--i >= 0) *pA++ = *pdata++;
				break;
							   }
			case mxINT8_CLASS:{
				char *pA=(char*)mxGetData(A);
				for(i=0;i<n;i++) *pA++ = data[i];
				//while (--i >= 0) *pA++ = *pdata++;
				break;
							  }
			default: assert(0);

		}
		return A;
	}


	void disp(){
		mxArray *A=tomxArray();
		MatlabInterface::disp(A);
		mxDestroyArray(A);
	}
	void ds(){
		mxArray *A=tomxArray();
		MatlabInterface::callMatlab("ds",A);
		mxDestroyArray(A);
	}
	void ds(int p,int q,int k){
		mxArray *A=tomxArray();
		setDimensions(A,p,q,k);
		MatlabInterface::callMatlab("ds",A);
		mxDestroyArray(A);
	}
	T sum(){
		T val=0;
		T*pdata=&data[0];
		int i=n;
		while (--i >= 0) val+=*pdata++;
		return val;
	}
	T max(){
		assert(n>0);
		return *max_element(data.begin(),data.end());
	}
	void transpose(){

		if(p>=2 && q>=2){
			vector<T>dataCopy(data);
			for(int j=0;j<q;j++){
				T*pdata=&data[j];
				T*pdataCopy=&dataCopy[j*p];
				for(int i=0;i<p;i++){
					*pdata=*pdataCopy++;
					pdata+=q;
				}
			}
		}
		/*		
		if(p>=2 && q>=2){
		T temp;
		for(int j=0;j<q;j++){
		T*pdata1=&data[j*p];
		T*pdata2=&data[j];
		for(int i=j+1;i<p;i++){
		temp=*pdata1;
		*pdata1++ = *pdata2;
		*pdata2=temp;
		pdata2+=p;
		}
		}
		}
		*/
		int pold=p;
		p=q;
		q=pold;
	}
	static void disp(vector<T>&v){
		Matrix<T>matrix;
		matrix.read_vector(v);
		matrix.transpose();
		mxArray *A=matrix.tomxArray();
		MatlabInterface::disp(A);
		mxDestroyArray(A);
	}
	static void deserialize_cell(const mxArray *A,vector<Matrix<T> >&v){
		int n=mxGetNumberOfElements(A);
		v.resize(n);
		for(int i=0;i<n;i++){
			v[i].readmxArray(mxGetCell(A,i));
		}
	}
	static mxArray* serialize_cell(vector<Matrix<T> >&v){
		return serialize_cell(v,MatlabInterface::getClassID(typeid(T)));
	}
	static mxArray* serialize_cell(vector<Matrix<T> >&v,mxClassID classID){
		int n=v.size();
		mxArray *A=mxCreateCellMatrix(n,1);		
		for(int i=0;i<n;i++){
			 mxSetCell(A,i,v[i].tomxArray(classID));
		}
		return A;
	}
	static void getField(const mxArray *A,int i,string field,Matrix<T>&matrix){
		matrix.readmxArray(getFieldByName(A,i,field));
	}
	static void getField(const mxArray *A,string field,Matrix<T>&matrix){
		getField(A,0,field,matrix);
	}
    static mxArray* createmxArray_empty(){
        int ndim=2;
		const mwSize dims[]={(mwSize)0,(mwSize)0};
        mxClassID classID=MatlabInterface::getClassID(typeid(T));
        mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
        return A;
    }
    static mxArray* createmxArray(int n){
        if(n==0)
            return createmxArray_empty();
        mxClassID classID=MatlabInterface::getClassID(typeid(T));
        int ndim=1;
		const mwSize dims[]={(mwSize)n};
        mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
        return A;
    }
    static mxArray* createmxArray(int p,int q){
        mxClassID classID=MatlabInterface::getClassID(typeid(T));
        int ndim=2;
		const mwSize dims[]={(mwSize)p,(mwSize)q};
        mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
        return A;
    }
    static mxArray* createmxArray(int p,int q,int k){
        mxClassID classID=MatlabInterface::getClassID(typeid(T));
        int ndim=3;
		const mwSize dims[]={(mwSize)p,(mwSize)q,(mwSize)k};
        mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
        return A;
    }
    
//     static mxArray* createmxArray_sparse(int p,int q,vector<int>&indi,vector<int>&indj,vector<T>&valij){
//         mxClassID classID=MatlabInterface::getClassID(typeid(T));
//         int ndim=3;
// 		const mwSize dims[]={(mwSize)p,(mwSize)q,(mwSize)k};
//         mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
//         return A;
//     }
    
    static mxArray* createmxArray_uninitialized(int n){
        if(n==0)
            return createmxArray_empty();        
        int ndim=1;
        const mwSize dims[]={(mwSize)0};
        mxClassID classID=MatlabInterface::getClassID(typeid(T));
        mxArray *A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
        void*data=mxGetData(A);
        mxFree(data);
        data=mxMalloc(n*sizeof(T));
        //             data=(double*)mxCalloc(n,sizeof(double));
        //             data=(double*)mxRealloc(data,n*sizeof(double));//no mxFree in that case
        //             double*pdata=data;
        //             int i=n;
        //             while (--i >= 0) *pdata++ = 0;
        mxSetData(A,data);
        
        int ndim2=1;
        const mwSize dims2[]={(mwSize)n};
        mxSetDimensions(A,dims2,ndim2);
        return A;        
    }
    static mxArray* createmxArray_uninitialized(int p,int q){
        mxArray*A=createmxArray_uninitialized(p*q);
        reshape(A,p,q);
        return A;        
    }
    static mxArray* createmxArray_uninitialized(int p,int q,int k){
        mxArray*A=createmxArray_uninitialized(p*q*k);
        reshape(A,p,q,k);
        return A;        
    }
    
    
    static void reshape(mxArray*A,int p,int q){
        int n=mxGetNumberOfElements(A);
        assert(p*q==n);
        int ndim=2;
		const mwSize dims[]={(mwSize)p,(mwSize)q};
         mxSetDimensions(A,dims,ndim);
    }
    static void reshape(mxArray*A,int p,int q,int k){
        int n=mxGetNumberOfElements(A);
        assert(p*q*k==n);
        int ndim=3;
		const mwSize dims[]={(mwSize)p,(mwSize)q,(mwSize)k};
         mxSetDimensions(A,dims,ndim);
    }
    
    
};

