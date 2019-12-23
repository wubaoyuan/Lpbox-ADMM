/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

mex_compatibility_mwSize
conversions for compatibility with mwSize
*=================================================================*/
#pragma once
//#include <assert.h>

#ifdef IS_NO_MWSIZE
typedef int mwSize;//unsigned int? size_t?
typedef int mwIndex;
#endif


inline int mwSize2int(mwSize m){
    int m2=(int)m;
    assert((mwSize)m2==m);//m<=pow(2,32)-1 ???
    return m2;
}
inline int mwIndex2int(mwIndex m){
    int m2=(int)m;
    assert((mwIndex)m2==m);//m<=pow(2,32)-1 ???
    return m2;
}
int mxGetM2(const mxArray *A){
    return mwSize2int(mxGetM(A));
}
int mxGetN2(const mxArray *A){
    return mwSize2int(mxGetN(A));
}
// const int*mxGetDimensions2(const mxArray *A){
//     const int *dims0 = mxGetDimensions(A);
//     int n=mxGetNumberOfDimensions(A);    
// }


