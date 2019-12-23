/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

mex_util = used by a couple of mex functions 
*=================================================================*/
//TODO:suppress disp()
//#undef NDEBUG

#pragma once
#include <sstream>
# include "math.h"
#include <string.h>
#include <time.h>
#include "mat.h"
#include <vector>
#include <map>
#include <string>
#include <cstdarg>
//#include <assert.h>



using namespace std;


//TODO:
// void compute_error_string(string condition,string line,string file,string &s){
//     s="MATLAB:error" + "!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )
//     s="In "+"<a href=""matlab: opentoline(''"+file+"''",',num2str(stack(i).line),');">',stack(i).name,' at ',num2str(stack(i).line),'</a>']);    
// }



//MACROS


//# define assert(isOK) ( (isOK) ? (void)0 : (void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )

#define printLine mexPrintf("line %d in file '%s'\n", __LINE__, __FILE__)
#define printError mexErrMsgIdAndTxt("MATLAB:error","line %d in file '%s'\n", __LINE__, __FILE__)

#ifndef NDEBUG
# define ASSERT( isOK ) \
	( (isOK) ? \
	(void)0 : \
	(void)mexErrMsgIdAndTxt("MATLAB:error","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )
//(void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )
#else
# define ASSERT( unused ) do {} while( false )
#endif

# define assert( isOK ) \
	( (isOK) ? \
	(void)0 : \
	(void)mexErrMsgIdAndTxt("MATLAB:error","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )
//(void)mexErrMsgIdAndTxt("","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) )

/*
# define assert2( isOK , msg) \
	( (isOK) ? \
	(void)0 : \
	(void) (mexErrMsgIdAndTxt("MATLAB:error","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__);disp2(msg,#msg) )
*/

/*
void assert(bool isOK ){
	if(!isOK)
		mexErrMsgIdAndTxt("MATLAB:error","!! Assert '%s' failed on line %d in file '%s'\n", #isOK, __LINE__, __FILE__) );
}
*/

// other inclusions
#include "MatlabInterface.cpp"
#include "mex_compatibility_mwSize.cpp"



//void disp(string s1,string s2);
void disp2(string s1,string s2);

#define DISP2(x)(disp2(x,#x))

#ifndef NDEBUG
# define DEBUG(val){mexPrintf("line %d, file '%s': ", __LINE__, __FILE__);DISP2(val);}
#else
//# define DEBUG(unused) {do {} while( false )}
# define DEBUG(unused) {}
#endif


void getSizes(const mxArray *A,int &p,int &q,int &k){
	if (mxGetNumberOfDimensions(A)<=2){
		p=mxGetM2(A);
		q=mxGetN2(A);
		k=1;
	}
	else {
		assert(mxGetNumberOfDimensions(A)<=3);
		const mwSize *dims = mxGetDimensions(A);
		p=mwSize2int(dims[0]);
		q=mwSize2int(dims[1]);
		k=mwSize2int(dims[2]);
	}
}
void getSizes(const mxArray *A,int &p,int &q,int &r,int&s){
	if(mxGetNumberOfDimensions(A)<=2){
		p=mxGetM2(A);
		q=mxGetN2(A);
		r=1;
        s=1;
	}
	else if(mxGetNumberOfDimensions(A)<=3){
		const mwSize *dims = mxGetDimensions(A);
		p=mwSize2int(dims[0]);
		q=mwSize2int(dims[1]);
		r=mwSize2int(dims[2]);
        s=1;
    }
    else{
		assert(mxGetNumberOfDimensions(A)<=4);
		const mwSize *dims = mxGetDimensions(A);
		p=mwSize2int(dims[0]);
		q=mwSize2int(dims[1]);
		r=mwSize2int(dims[2]);
		s=mwSize2int(dims[3]);
	}
}
int getSize3(const mxArray *A){
	int k=0;
	if (mxGetNumberOfDimensions(A)<=2)
		k = 1; 
	else {
		const mwSize *dims = mxGetDimensions(A);
		k = mwSize2int(dims[2]);
	}
	return k;
}

void** deal_output(int num,...) {
	//deprecated:use deal_outputs
	void**output=(void**)mxCalloc(num,sizeof(void*));
	va_list arguments;                     
	va_start(arguments,num);           
	for(int i=0;i<num;i++)        
		output[i]=va_arg(arguments,void*); 
	va_end(arguments);                  
	return output;
}
vector<void*> deal_outputs(int num,...) {
	vector<void*>output(num);
	//void**output=(void**)mxCalloc(num,sizeof(void*));
	va_list arguments;                     
	va_start(arguments,num);           
	for(int i=0;i<num;i++)        
		output[i]=va_arg(arguments,void*); 
	va_end(arguments);                  
	return output;
}
mxArray* createMxArrayInt(int m,int n){
	int ndim=2;
	const mwSize dims[]={(mwSize)m,(mwSize)n};
	return mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
}
mxArray* createMxArrayInt(int m,int n,int k){
	int ndim=3;
	const mwSize dims[]={(mwSize)m,(mwSize)n,(mwSize)k};
	return mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
}

void fill_out_full_int(mxArray *out[],int nout,int *val,int m,int n) {
	out[nout] = mxCreateDoubleMatrix(m, n, mxREAL);
	double *temp = mxGetPr(out[nout]);
	int nnz = m*n;
	for (int i=0;i<nnz;i++)
		temp[i] = (double)val[i];
}
void fill_out_full(mxArray *out[],int nout,double *val,int m,int n) {
	out[nout] = mxCreateDoubleMatrix(m, n, mxREAL);
	double *temp = mxGetPr(out[nout]);
	int nnz = m*n;
	for (int i=0;i<nnz;i++)
		temp[i] = val[i];
}

void readSparse(const mxArray *A,double**pr,mwIndex**ir,mwIndex**jc,mwSize *m,mwSize *n){
    assert(mxIsDouble(A));
	*pr=mxGetPr(A);
	*ir=mxGetIr(A);
	*jc=mxGetJc(A);
	*m=mxGetM(A);
	*n=mxGetN(A);
}
void readSparse(const mxArray *A,bool**pr,mwIndex**ir,mwIndex**jc,mwSize *m,mwSize *n){
    assert(mxIsLogical(A));
	*pr=(bool*)mxGetData(A);
	*ir=mxGetIr(A);
	*jc=mxGetJc(A);
	*m=mxGetM(A);
	*n=mxGetN(A);
}

int* readInt32(const mxArray *A) {
	switch (mxGetClassID(A))  {
		case mxINT32_CLASS: return (int*)mxGetData(A);  break;//TODO:attention with subsequent modif of pointer
		case mxDOUBLE_CLASS: {
			double *pr = mxGetPr(A);
			int n = mxGetM2(A)*mxGetN2(A);
			int *res;
			if(n>0){
				//assert(0); memory leak ?
				res = (int*)mxCalloc(n,sizeof(int));
			}
			else
				res = NULL;
			for (int i=0;i<n;i++)
				res[i] = (int)pr[i];
			return res;
							 }
		default: assert(0);//mexErrMsgTxt("bad input");
			return NULL;
	}
}

vector<string>getFieldNames(const mxArray*A){
	int n=mxGetNumberOfFields(A);
	vector<string> fields;
	fields.resize(n);
	for(int i=0;i<n;i++)
		fields[i]=mxGetFieldNameByNumber(A,i);
	return fields;
}
bool isField(const mxArray*A,string field){
	vector<string>fields=getFieldNames(A);
	for(int i=0;i<fields.size();i++)
		if(fields[i]==field)
			return true;
	return false;
}
//mxArray*getFieldByName(const mxArray*A,int i, const char *fname){
mxArray*getFieldByName(const mxArray*A,int i, string fname){
	assert(mxIsStruct(A)&&isField(A,fname));//inefficient?
	int field=mxGetFieldNumber(A,fname.c_str());
	return mxGetFieldByNumber(A,i,field);
}
mxArray*getFieldByName(const mxArray*A,string fname){
	return getFieldByName(A,0, fname);
}
//void setFieldByName(mxArray*A,int i, const char *fname,mxArray*B){
void setFieldByName(mxArray*A,int i, string fname,mxArray*B){
	assert(isField(A,fname));//inefficient?
	int field=mxGetFieldNumber(A,fname.c_str());
	mxSetFieldByNumber(A,i,field,B);      
}
mxArray* array2mxArray(double*x,int m,int n) {
	mxArray* A=mxCreateDoubleMatrix(m, n, mxREAL);
	double *pA=mxGetPr(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* array2mxArray(bool*x,int m,int n) {//TODO:make better (but obsolete anyway)
	mxArray* A=mxCreateLogicalMatrix(m, n);
	bool *pA=(bool*)mxGetData(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* array2mxArray(int*x,int m,int n) {//TODO:make better (but obsolete anyway)
	int ndim=2;
	const mwSize dims[]={(mwSize)m,(mwSize)n};
	mxArray* A=mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL);
	int *pA=(int*)mxGetData(A);
	int mn=m*n;    
	for(int i=0;i<mn;i++,pA++,x++)
		*pA=*x;
	x-=mn;
	return A;
}
mxArray* array2mxArray(int*x,int m,int n,mxClassID classID) {//TODO:make better (but obsolete anyway)
	int ndim=2;
	const mwSize dims[]={(mwSize)m,(mwSize)n};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	int mn=m*n;    
	switch(classID){
		case mxDOUBLE_CLASS:{
			double *pA=mxGetPr(A);
			for(int i=0;i<mn;i++,pA++,x++)
				*pA=(double)*x;
			x-=mn;
			break;
							}
		default:
			assert(0);//mexErrMsgTxt("wrong argument\n");
	}
	return A;
}
mxArray* array2mxArray(double*x,int m,int n,int k,mxClassID classID) {
	int ndim=3;
	const mwSize dims[]={(mwSize)m,(mwSize)n,(mwSize)k};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	int mnk=m*n*k;    
	switch(classID){
		case mxDOUBLE_CLASS:{
			double *pA=mxGetPr(A);
			for(int i=0;i<mnk;i++,pA++,x++)
				*pA=(double)*x;
			x-=mnk;
			break;
							}
		default:
			assert(0);//mexErrMsgTxt("wrong argument\n");
	}
	return A;
}


mxArray* array2mxArray(vector<double> &v){
	int n=v.size();
	mxArray* A;
	//if(n==0)
	//	A=mxCreateDoubleMatrix(0, 0, mxREAL);
	//else
	//	A=mxCreateDoubleMatrix(n, 1, mxREAL);

	A=mxCreateDoubleMatrix(n, 1, mxREAL);
	double *pA=mxGetPr(A);
	for(int i=0;i<n;i++)
		*pA++=v[i];
	return A;
}
mxArray* array2mxArray(vector<bool> &v){
	int n=v.size();
	mxArray* A=mxCreateLogicalMatrix(n, 1);
	bool *pA=(bool*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pA++=v[i];
	return A;
}
mxArray* array2mxArray(vector<int> &v){
	int n=v.size();
	mxArray* A=createMxArrayInt(n, 1);
	int *pA=(int*)mxGetData(A);
	int *pv=&v[0];
	for(int i=0;i<n;i++)
		*pA++=*pv++;
	return A;
}
mxArray* array2mxArray(vector<int> &v,mxClassID classID,int a){
	//mxArray*A=array2mxArray(v,mxDOUBLE_CLASS,0);
	//mxDOUBLE_CLASS,mxINT32_CLASS
	int n=v.size();
	vector<int>::iterator pv=v.begin();
	switch (classID)  {
		case mxINT32_CLASS:			
			{
				mxArray* A=createMxArrayInt(n,1);
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pA++=a+*pv++;
				return A;  
				break;
			}						
		case mxDOUBLE_CLASS: 
			{
				double a2=(double)a;
				mxArray* A=mxCreateDoubleMatrix(n,1,mxREAL);
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pA++=a2+(double)*pv++;
				return A;  
				break;										
			}			
		default: assert(0);//mexErrMsgTxt("bad input");
			return NULL;
	}
}
mxArray* array2mxArray(vector<char> &v){
	int n=v.size();
	vector<char>::iterator pv=v.begin();
	int ndim=2;
	const mwSize dims[]={(mwSize)n,(mwSize)1};
	mxArray* A=mxCreateNumericArray(ndim,dims,mxINT8_CLASS,mxREAL);

	char*pA=(char*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pA++=*pv++;
	return A;  
}
mxArray* array2mxArray(vector<unsigned char> &v){
	int n=v.size();
	vector<unsigned char>::iterator pv=v.begin();
	int ndim=2;
	const mwSize dims[]={(mwSize)n,(mwSize)1};
	mxArray* A=mxCreateNumericArray(ndim,dims,mxUINT8_CLASS,mxREAL);

	unsigned char*pA=(unsigned char*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pA++=*pv++;
	return A;  
}
void mxArray2array(const mxArray* A,vector<char> &v){
	int n = mxGetNumberOfElements(A);
	v.resize(n);
	vector<char>::iterator pv = v.begin();
	assert(mxGetClassID(A) == mxINT8_CLASS);
	char* pA = (char*) mxGetData(A);
	for(int i = 0; i < n; i++)
		*pv++=*pA++;			
}
void mxArray2array(const mxArray* A,vector<unsigned char> &v){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<unsigned char>::iterator pv=v.begin();
	assert(mxGetClassID(A)==mxUINT8_CLASS);
	unsigned char*pA=(unsigned char*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pv++=*pA++;			
}
void mxArray2array(const mxArray* A,vector<int> &v,int a){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<int>::iterator pv=v.begin();
	switch (mxGetClassID(A)) {
		case mxINT32_CLASS: 
			{
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pv++=(a+*pA++);			
				break;							
			}
		case mxDOUBLE_CLASS: 
			{
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pv++=(int)(a+*pA++);
				break;							
			}

		default: assert(0);//mexErrMsgTxt("bad input");
	}	
}
void mxArray2array(const mxArray* A,vector<double> &v,double a){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	vector<double>::iterator pv=v.begin();
	switch (mxGetClassID(A)) {
		case mxINT32_CLASS: 
			{
				int*pA=(int*)mxGetData(A);
				for(int i=0;i<n;i++)
					*pv++=(double)(a+*pA++);			
				break;							
			}
		case mxDOUBLE_CLASS: 
			{
				double*pA=mxGetPr(A);
				for(int i=0;i<n;i++)
					*pv++=(double)(a+*pA++);
				break;							
			}
		default: printError;
	}	
}
void mxArray2array(const mxArray* A,vector<bool> &v){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	//vector<bool>::iterator pv=v.begin();
	assert(mxIsLogical(A));
	bool*pA=(bool*)mxGetData(A);
	for(int i=0;i<n;i++)
		v[i]=*pA++;	
		//*pv++=*pA++;	
	//switch (mxGetClassID(A)) {
	//	case mxLOGICAL_CLASS:
	//		{
	//			bool*pA=(bool*)mxGetData(A);
	//			for(int i=0;i<n;i++)
	//				*pv++=*pA++;	
	//			break;											
	//		}
	//	default: printError;
	//}	
}
void mxArray2array(const mxArray* A,vector<int> &v){
	int n=mxGetNumberOfElements(A);
	v.resize(n);
	int*pv=&v[0];
	assert(mxGetClassID(A)==mxINT32_CLASS);
	int*pA=(int*)mxGetData(A);
	for(int i=0;i<n;i++)
		*pv++=*pA++;	
}
void mxArray2array(const mxArray* A,string &v){
	assert(mxGetClassID(A)==mxCHAR_CLASS);
	char *v0=mxArrayToString(A);
	v=v0;
}
void mxArrayStruct2array(const mxArray* A,int i,string field,string &v){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,vector<double> &v,double a){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v,a);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,vector<int> &v,int a){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v,a);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,vector<bool> &v){
	const mxArray*Ai=getFieldByName(A,i,field);
	mxArray2array(Ai,v);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,int &v,int a){
	const mxArray*Ai=getFieldByName(A,i,field);
	v=a+(int)mxGetScalar(Ai);
}
void mxArrayStruct2array(const mxArray* A,int i,string field,double &v,double a){
	const mxArray*Ai=getFieldByName(A,i,field);
	v=a+(double)mxGetScalar(Ai);
}


template<class T>int sizeVector(vector<T> &v){
	return v.size()*sizeof(T);
}
mxArray* create_mxArray(int m,int n,mxClassID classID) {
	int ndim=2;
	const mwSize dims[]={(mwSize)m,(mwSize)n};
	mxArray* A=mxCreateNumericArray(ndim,dims,classID,mxREAL);
	return A;
}
double* createArray(int n,double val){
	if(n==0)
		return NULL;
	double *x = (double*)mxCalloc(n,sizeof(double));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}

bool* createArray(int n,bool val){
	if(n==0)
		return NULL;
	bool *x = (bool*)mxCalloc(n,sizeof(bool));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}
int* createArray(int n,int val){
	if(n==0)
		return NULL;
	int *x = (int*)mxCalloc(n,sizeof(int));
	for(int i=0;i!=n;i++,x++)
		*x=val;
	x-=n;
	return x;
}

mxArray*createStructArray(int n, vector<string>fields){
	int nFields=fields.size();
	const char **field_names = (const char **)mxCalloc(nFields, sizeof(*field_names));
	for(int i=0;i<nFields;i++)
		field_names[i]=fields[i].c_str();
	mxArray*A=mxCreateStructMatrix(n,1,nFields,field_names);
	mxFree(field_names);
	return A;
}
//bool isEmpty(const mxArray*A){
//	return mxGetNumberOfElements(A)==0;
//}


void compute_cell_strings(const mxArray *A,vector<string>&strings,bool isValid_single){
    if(mxIsCell(A)){
        int n=mxGetNumberOfElements(A);
        strings.resize(n);        
        for(int s=0;s<n;s++){
            const mxArray *As=mxGetCell(A,s);
            if(As==NULL)
                assert(0);
            assert(mxGetClassID(As)==mxCHAR_CLASS);
            strings[s]=mxArrayToString(As);
        }
    }
    else{
        assert(isValid_single);
        if(A==NULL)       
            strings.resize(0);
        else{
            assert(mxGetClassID(A)==mxCHAR_CLASS);
            strings.resize(1);
            strings[0]=mxArrayToString(A);
        }
    }
}

/* use mxArrayToString
char *str readString(const mxArray *A) {
int buflen = (mxGetM2(A) * mxGetN2(A)) + 1;
char *str=mxCalloc(buflen, sizeof(char));
int status = mxGetString(A, str, buflen);
return str;
}
*/

void copyDims(const mxArray *A,mxArray *B){
    assert(mxGetNumberOfElements(A)==mxGetNumberOfElements(B));//TODO:voir if we need that
	mxSetDimensions(B,mxGetDimensions(A), mxGetNumberOfDimensions(A));
}
void setDimensions(mxArray *A,int m,int n){
	int ndim=2;
	const mwSize dims[]={(mwSize)m,(mwSize)n};
	assert(mxGetNumberOfElements(A)==m*n);
	int isError = mxSetDimensions(A, dims, ndim);
	if(isError)
		printError;
}

void setDimensions(mxArray *A,int m,int n,int k){
	int ndim=3;
	const mwSize dims[]={(mwSize)m,(mwSize)n,(mwSize)k};
	assert(mxGetNumberOfElements(A)==m*n*k);
	int isError = mxSetDimensions(A, dims, ndim);
	if(isError)
		printError;
}

void swap(double *x,double *y) {
	double temp=*x;
	*x=*y;
	*y=temp;
}
void swap(double&x,double&y) {
	double temp=x;
	x=y;
	y=temp;
}

double *copyArray(const mxArray *A) {
	//TODO:use mxDuplicateArray?
	int n=mxGetM2(A)*mxGetN2(A);
	return (double*) memcpy(mxCalloc(n, sizeof(double)), mxGetPr(A), sizeof(double)*n);
}
double *copyArray(double *x,int n) {
	return (double*) memcpy(mxCalloc(n, sizeof(double)), x, sizeof(double)*n);
}
void copyArray(double *x,int n,double *y) {
	memcpy(y, x, sizeof(double)*n);
}


void compute_j2Jc(mwIndex*jc,int*indj,int n,int nnz) {    
	//assumes indi,indj ordered
	int k = 0;    
	for(int j=0;j<n;j++) {
		jc[j] = (mwIndex)k;
		while(k<nnz && indj[k]==j)
			k++;
	}
	jc[n] = (mwIndex)k;
}
void fill_irjc_from_ind(int*ind,int*ir,int*jc,int nnz,int p,int q){
	//assumes ind sorted
	int k = 0;    
	int u,ju;
	//int*pind=ind;
	for(int j=0;j<q;j++) {
		jc[j] = k;
		while(k<nnz){
			u=*ind;
			ju=u/p;
			*ir=u-ju*p;
			if(ju==j){
				k++;
				ind++;
				ir++;
			}
			else
				break;
		}
		//while(k<nnz && ind[k]==j)
		//	k++;
	}
	jc[q] = k;
	//for(int u=0;u<nnz;u++)
	//	*ir++=*ind++%p;
}
/*
int compute_ij2k(int*jc,int*indj,int n,int nnz) {
int k = 0;    
for(int j=0;j<n;j++) {
jc[j] = k;
while(k<nnz && indj[k]==j) {
k++;
}        
}
jc[n] = k;
}
*/

mxArray* sparseSorted(int*indi,int*indj,double*val,int m,int n,int nnz){
	mxArray*W=mxCreateSparse(m,n,nnz,mxREAL);
	double *pr = mxGetPr(W);
	mwIndex *ir = mxGetIr(W);
	mwIndex *jc = mxGetJc(W);
	int k;
	mwIndex*pir=ir;
    int*pindi=indi;
	double*ppr=pr,*pval=val;
	for(k=0;k<nnz;k++){
		*pir++=(mwIndex)*pindi++;
		*ppr++=*pval++;
	}
	compute_j2Jc(jc,indj,n,nnz);
	return W;
}
mxArray* sparseSorted(int*indi,int*indj,double val,int m,int n,int nnz){
	mxArray*W=mxCreateSparse(m,n,nnz,mxREAL);
	double *pr = mxGetPr(W);
	mwIndex *ir = mxGetIr(W);
	mwIndex *jc = mxGetJc(W);
	int k;
	mwIndex*pir=ir;
    int*pindi=indi;
	double*ppr=pr;
	for(k=0;k<nnz;k++){
		*pir++=(mwIndex)*pindi++;
		*ppr++=val;
	}
	compute_j2Jc(jc,indj,n,nnz);
	return W;
}

/*
void fill_out_sparse(mxArray *out[],int nout,int *pr, int *jc,double *pr,int nnz, int m,int n) {
out[nout] = mxCreateSparse(m, n, nnz, mxREAL);
double *temp = mxGetPr(out[nout]);
for (int i=0;i<nnz;i++)
temp[i] = val[i];
}*/


void sparse2full(const mxArray* A,vector<double>&x){
	assert(mxIsSparse(A));
	int m=mxGetM2(A);
	int n=mxGetN2(A);
	int N=m*n;
	x.clear();//remove ?
	x.resize(N,0);
	mwIndex*ir=mxGetIr(A);
	mwIndex*jc=mxGetJc(A);
	double*pr=mxGetPr(A);
	int k,i;
	for (int j=0;j<n;j++)
		for (k=(mwIndex)jc[j]; k!=(mwIndex)jc[j+1]; k++) {
			i = (mwIndex)ir[k];
			x[i+j*m]=pr[k];
		}
}
void displayMatrix(double *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = W[k];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 1, args, "disp");
	mxDestroyArray(W2);
}
void displayMatrix(int *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = (double)(W[k]);
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	//NOTE: syntax will not work with -argcheck compiler option
	mexCallMATLAB(0, NULL, 1, args, "disp");
	mxDestroyArray(W2);
}

void imagescMatrix(double *W,int m,int n) {
	mxArray *W2 = mxCreateDoubleMatrix(m, n, mxREAL);
	double *prW = mxGetPr(W2);
	int k;
	for (k=0;k<m*n;k++)
		prW[k] = W[k];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 0, NULL, "figure2");
	mexCallMATLAB(0, NULL, 1, args, "imagesc");
	mxDestroyArray(W2);
}
void showImages(double *W,int m,int n,int k) {
	const mwSize dims[]={(mwSize)m,(mwSize)n,(mwSize)k};
	mxArray *W2 = mxCreateNumericArray(3,dims,  mxDOUBLE_CLASS,mxREAL);
	double *prW = mxGetPr(W2);
	int i;
	for (i=0;i<m*n*k;i++)
		prW[i] = W[i];
	mxArray *args[1];
	args[0] = (mxArray*) W2;
	mexCallMATLAB(0, NULL, 0, NULL, "figure2");
	mexCallMATLAB(0, NULL, 1, args, "showImages");
	mxDestroyArray(W2);
}


void debug(bool condition,int val) {
	if(condition) {
		mexPrintf("val1 = %d\n",val);
	}
}
void debug(bool condition,int val,char *s) {
	if(condition) {
		//mexPrintf(s+" = %d\n",val);
		//mexPrintf(strcat(s," = %d\n"),val);
		printf(strcat(s," = %d\n"),val);
	}
}
void debug(bool condition,double val) {
	if(condition) {
		mexPrintf("val1 = %1.3g\n",val);
	}
}
void stop() {
	//mexErrMsgTxt("Stop\n");
	assert(0);
}
void stop(char *s) {
	mexPrintf("%s\n",s);    
	//mexErrMsgTxt("Stop\n");
	assert(0);
}
//void stop(int i,char *s) {
//	//stop(__LINE__,__FILE__);
//	//mexPrintf("Stopped at: File: %s, line: %d\n",__FILE__,__LINE__);
//	mexPrintf("Stopped at: File: %s, line: %d\n",s,i);
//	mexErrMsgTxt("Stop\n");
//}

void checkBounds(int val,int lim) {
	if(val<0) {		
		mexPrintf("checkBounds failed: val = %d < 0\n",val);
		assert(0);//mexErrMsgTxt("Error\n");
	}
	if(val>=lim){
		mexPrintf("checkBounds failed: val = %d >= %d\n",val,lim);
		assert(0);//mexErrMsgTxt("Error\n");
	}
}

void disp(string s,string var) {
	mexPrintf("%s=%s\n",var.c_str(),s.c_str());
}
void disp(string s) {
	mexPrintf("%s\n",s.c_str());
}
void disp(int x) {
	mexPrintf("? = %d\n",x);
}
void disp(double x) {
	//mexPrintf("? = %1.3g\n",x);
	mexPrintf("? = %1.6g\n",x);
}
void disp(int x,char *s) {
	mexPrintf("%s = %d\n",s,x);
	//mexPrintf(strcat(s," = %d\n"),x);
}
void disp(double x,char *s) {
	//mexPrintf("%s = %1.3g\n",s,x);
	mexPrintf("%s = %1.6g\n",s,x);
	//mexPrintf(strcat(s," = %1.3g\n"),x);
}


void disp(vector<double> &x){
	disp("vector:");
	for(int i=0;i<x.size();i++)
		mexPrintf("?[%d] = %1.3g\n",i,x[i]);
}
void disp(vector<int> &x){
	disp("vector:");
	for(int i=0;i<x.size();i++)
		mexPrintf("?[%d] = %d\n",i,x[i]);
}
void disp(vector<vector<int> > &x){
	disp("vector:");
	for(int j=0;j<x.size();j++){
		mexPrintf("?[%d]:",j);
		for(int i=0;i<x[j].size();i++)
			mexPrintf("?[%d][%d] = %d ; ",j,i,x[j][i]);
		mexPrintf("\n");
	}
}
template<class T>void disp2(T val,string var){
	const type_info& t = typeid(T);	
	if(t==typeid(int) || t==typeid(char) || t==typeid(unsigned int) || t==typeid(unsigned char)  || t==typeid(typename vector<T>::size_type) || t==typeid(bool)){
		mexPrintf("%s = %d\n",var.c_str(),(int)val);		
	}
	else if(t==typeid(double))
		mexPrintf("%s = %1.6g\n",var.c_str(),val);
	else if(t==typeid(float))
		mexPrintf("%s = %1.6g\n",var.c_str(),val);
	//else if(t==typeid(bool))
	//	mexPrintf("%s = %d\n",var.c_str(),(int)val);
	else if(t==typeid(char *))
		//mexPrintf("%s = %s\n",var.c_str(),((string)val).c_str());
		mexPrintf("%s = %s\n",var.c_str(),val);
// 	else if(t==typeid(string))
// 		mexPrintf("%s = %s\n",var.c_str(),val.c_str());//compile error/warning
	else {
		assert(0);
	}
}
void disp2(string val,string var){
	mexPrintf("%s = %s\n",var.c_str(),val.c_str());
}

template<class T>void disp2(T val){
	disp2(val,"?");
}

void dispTime(long time){
	mexPrintf("time = %10.5g\n",((double)time)/CLOCKS_PER_SEC);   
}
void dispTime(long time[],int n){
	double sum=0,timei;
	for(int i=0;i<n;i++){
		timei=((double)time[i])/CLOCKS_PER_SEC;
		mexPrintf("time[%d] = %10.5g\n",i,timei);   
		sum+=timei;
	}
	mexPrintf("TOTAL = %10.5g\n",sum);   
}
void dispTime(map<string,long>&timing){
	double sum=0,timei;

	for(map<string,long>::iterator iter = timing.begin(); iter != timing.end(); iter++){
		timei=((double)(*iter).second)/CLOCKS_PER_SEC;
		mexPrintf("timing[%s] = %10.5g\n",(*iter).first.c_str(),timei);   
		sum+=timei;
	}
	mexPrintf("TOTAL = %10.5g\n",sum);   
}
template<class T>void deleteVector(vector<T> &v){
	//TODO:v.clear(); as well ?
	int n=v.size();
	for(int i=0;i<n;i++){
		assert(v[i]!=NULL);
		delete v[i];
	}
}
template<class T>void array2vector(T*x,vector<T> &v,int n){
	v.resize(n);
	typename vector<T>::iterator pv=v.begin();
	for(int i=0;i<n;i++)
		*pv++=*x++;
}
template<class T>void vector2array(vector<T> &v,T*x){
	int n=v.size();
	typename vector<T>::iterator pv=v.begin();
	for(int i=0;i<n;i++)
		*x++=*pv++;
}





template<class T>void disp_vector(vector<T> &x){
	disp("vector:");
	const type_info& t = typeid(T);	
	if(t==typeid(int) || t==typeid(char) || t==typeid(unsigned int) || t==typeid(unsigned char)  || t==typeid(bool) || t==typeid(typename vector<T>::size_type) )
		for(int i=0;i<x.size();i++)
			mexPrintf("?[%d] = %d\n",i,(int)x[i]);
	else if(t==typeid(double))
		for(int i=0;i<x.size();i++)
			mexPrintf("?[%d] = %1.6g\n",i,x[i]);
	else {
		assert(0);
	}
}
string toStr(int x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
string toStr(double x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
string toStr(long x){
	ostringstream s0;
	s0 << x;
	string s = s0.str();
	return s;
}
/*
string formatFilename(string file){
string file2;
for(int i=0;i<file.size();i++){
if(file[i]=='\')
file2.push_back('\\');
file2.push_back(file[i]);
}
}
*/



string filepath_save_trash="C:\\tim\\temp\\trash\\";

void TEMPSAVE(mxArray*mxA,string name,int p,int q,int k){
	setDimensions(mxA,p,q,k);
	MatlabInterface::save2matfile(mxA,filepath_save_trash+name+".mat",name);
}
void TEMPSAVE(mxArray*mxA,string name,int p,int q){
	setDimensions(mxA,p,q);
	MatlabInterface::save2matfile(mxA,filepath_save_trash+name+".mat",name);
}
void TEMPSAVE(vector<int>&v,string name,int p,int q,int k){
	mxArray*mxA=array2mxArray(v,mxDOUBLE_CLASS,0);
	setDimensions(mxA,p,q,k);
	MatlabInterface::save2matfile(mxA,filepath_save_trash+name+".mat",name);
}
void TEMPSAVE(vector<int>&v,string name,int p,int q){
	//TEMPSAVE(v,"A",p,q);
	mxArray*mxA=array2mxArray(v,mxDOUBLE_CLASS,0);
	setDimensions(mxA,p,q);
	MatlabInterface::save2matfile(mxA,filepath_save_trash+name+".mat",name);
}
void TEMPSAVE(mxArray*X,string name){
	MatlabInterface::save2matfile(X,filepath_save_trash+name+".mat",name);
}
