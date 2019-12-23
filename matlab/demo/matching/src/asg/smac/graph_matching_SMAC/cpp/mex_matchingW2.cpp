//function [W,indi1i2j1j2,indiW_indjW] = mex_matchingW2(G1,W1,G2,W2,F12,I12,W12,f,options);
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.


# undef NDEBUG
# include "GraphMatching.cpp"
void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
  int k = 0;
  const mxArray *mxArray_G1 = in[k++];
  const mxArray *mxArray_W1 = in[k++];
  const mxArray *mxArray_G2 = in[k++];
  const mxArray *mxArray_W2 = in[k++];
  const mxArray *mxArray_F12 = in[k++];
  const mxArray *mxArray_I12 = in[k++];
  const mxArray *mxArray_W12 = in[k++];
  const mxArray *mxArray_f = in[k++];
  const mxArray *options = in[k++];

  GraphMatching *graph = new GraphMatching();
  graph->deserializeGraphs(mxArray_W1, mxArray_W2, mxArray_I12, mxArray_W12, mxArray_G1, mxArray_G2, mxArray_F12, mxArray_f);
  graph->deserializeOptions(options);
  graph->compute_W();
  graph->normalize_W();

  out[0] = graph->serialize_W();
  if(nargout >= 2)
    out[1] = graph->serialize_indi1i2j1j2();
  if(nargout >= 3)
    out[2] = graph->serialize_indiW_indjW();
  
  delete graph;
}

