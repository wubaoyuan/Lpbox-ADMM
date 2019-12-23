//function [indexes,L] = mex_classes2indexes(classes,k);
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.


#undef NDEBUG
#include <mex.h>
#include "mex_util.cpp"
#include <algorithm>
#include "Matrix.cpp"

void classes2indexes(Matrix<int>& classes,vector<Matrix<int> >& indexes,Matrix<int>& L){
	int k=L.p;
	int i,j;
	int*pclasses=&classes.data[0];
	for(i=0;i<classes.n;i++){
		if((j=pclasses[i])>=0)
			L.data[j]++;
	}
	for(j=0;j<k;j++){
		indexes[j].resize(L.data[j],1);
	}
	vector<int>temp(k);
	pclasses=&classes.data[0];
	for(i=0;i<classes.n;i++){
		//j=*pclasses++;
		j=pclasses[i];
		//indexes[j].data[temp[j]++]=i;
		//indexes[j].data[temp[j]++]=i;
		if(j>=0)
			indexes[j][temp[j]++]=i;
	}
}
void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
	Matrix<int> classes(in[0]);
	int k;
	if(nargin<2){
		if(classes.n>=0)
			k=classes.max();
		else
			k=0;
	}
	else
		k=(int)mxGetScalar(in[1]);
	assert(k>=0);//0 meaning indexes={};
	classes.add(-1);
	Matrix<int> L(k,1);
	vector<Matrix<int> > indexes(k);
	classes2indexes(classes,indexes,L);
	for(int i=0;i<k;i++){
		indexes[i].add(1);
	}
	out[0]=Matrix<int>::serialize_cell(indexes,mxDOUBLE_CLASS);
	out[1]=L.tomxArray(mxDOUBLE_CLASS);
}


