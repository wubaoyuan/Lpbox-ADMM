// function k2=mex_find_next_nonzero(L,ind,k);
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

#undef NDEBUG

#include "math.h"
#include "mex.h"
#include <vector>
#include "mex_util.cpp"

int find_next_nonzero(bool* L, int* ind, int k){
    //int n=ind.size();
    int* pind = ind + k;
    bool* pL = L - 1;
    while (!*(pL + *pind++)) {
	k++;	
	//assert(k<n);
    }
    //while(!L[*pind++]){
    //	k++;	
    //	//assert(k<n);
    //}
    return k;
}

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    //vector<bool>L;
    //vector<int>ind;
    //mxArray2array(in[0],L);
    //mxArray2array(in[1],ind,-1);

    bool* L = (bool*) mxGetData(in[0]);
    int* ind = (int*) mxGetData(in[1]);

    int k = (int) mxGetScalar(in[2]) - 1;
    int k2 = find_next_nonzero(L, ind, k);
    out[0] = mxCreateDoubleScalar((double) k2 + 1);
}
