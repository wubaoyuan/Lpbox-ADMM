// class MatlabInterface
// Timothee Cour, 04-Aug-2008 20:46:38 -- DO NOT DISTRIBUTE


#pragma once
#include <typeinfo>

class MatlabInterface {
public: 
	static mxClassID getClassID(const type_info& t){
		//mxClassID classID=getClassID(typeid(T));
		//const type_info& t = typeid(T);
		if(t==typeid(double)) return mxDOUBLE_CLASS;
		if(t==typeid(float)) return mxSINGLE_CLASS;
		if(t==typeid(int)) return mxINT32_CLASS;
		if(t==typeid(unsigned int)) return mxUINT32_CLASS;
		if(t==typeid(unsigned char)) return mxUINT8_CLASS;
		if(t==typeid(char)) return mxINT8_CLASS;
		if(t==typeid(bool)) return mxLOGICAL_CLASS;
		assert(0);
		return mxUNKNOWN_CLASS;
	}
	static void disp(mxArray*X){
		callMatlab("disp",X);
	}

	static void callMatlab(mxArray*&Y1,string command,mxArray*X1){
		const int nlhs=1;
		const int nrhs=1;
		mxArray *plhs[nlhs];
		mxArray *prhs[nrhs];
		prhs[0]=X1;
		int res=mexCallMATLAB(nlhs,plhs,nrhs,prhs,command.c_str());
		assert(res==0);
		Y1=plhs[0];
	}
	static void callMatlab(mxArray*&Y1,string command,mxArray*X1,mxArray*X2){
		const int nlhs=1;
		const int nrhs=2;
		mxArray *plhs[nlhs];
		mxArray *prhs[nrhs];
		prhs[0]=X1;
		prhs[1]=X2;
		int res=mexCallMATLAB(nlhs,plhs,nrhs,prhs,command.c_str());
		assert(res==0);
		Y1=plhs[0];
	}
	static void callMatlab(string command,mxArray*X1){
		const int nlhs=0;
		const int nrhs=1;
		mxArray *plhs[1];//?
		mxArray *prhs[nrhs];
		prhs[0]=X1;
		int res=mexCallMATLAB(nlhs,plhs,nrhs,prhs,command.c_str());
		assert(res==0);
	}
	static void callMatlab(string command){
		const int nlhs=0;
		const int nrhs=0;
		mxArray *plhs[1];//?
		mxArray *prhs[1];//?
		int res=mexCallMATLAB(nlhs,plhs,nrhs,prhs,command.c_str());
		assert(res==0);
	}
	static void save2matfile(mxArray *A,string file,string varname){
		//mxArray*mxA=array2mxArray(scores);
		//setDimensions(mxA,m,n);
		//save2matfile(mxA,"C:\\tim\\temp\\trash\\1.mat","A");
        
        //matOpen doesn't work under windows unless some library is included
        #if defined(_WIN32)
        assert(0);
        #else
        MATFile *pmat = matOpen(file.c_str(), "w");
        if (pmat == NULL)
            assert(0);//mexErrMsgTxt("error accessing mat file\n");
        int status = matPutVariable(pmat, varname.c_str(), A);
        if (status != 0)
            assert(0);//mexErrMsgTxt("error writing mat file\n");
        if (matClose(pmat) != 0)
            mexPrintf("Error closing file %s\n",file.c_str());
        else
            mexPrintf("saved file %s\n",file.c_str());
        #endif
	}
};

