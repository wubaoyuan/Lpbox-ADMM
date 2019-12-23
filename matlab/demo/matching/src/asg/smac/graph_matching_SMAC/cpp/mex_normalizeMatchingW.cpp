// function [W, D1, D2] = mex_normalizeMatchingW(W, I12, options);

# include "GraphMatching.cpp"
#undef NDEBUG

void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    int k = 0;
    const mxArray *mxArray_W = in[k++];
    const mxArray *mxArray_I12 = in[k++];
    const mxArray *options = in[k++];

    GraphMatching *graphMatching = new GraphMatching();
    graphMatching->deserializeOptions(options);
    graphMatching->deserializeW(mxArray_W, mxArray_I12);
    graphMatching->normalize_W();

    out[0] = graphMatching->serialize_W();
    if (nargout >= 2)
	out[1] = graphMatching->serialize_D1();
    if (nargout >= 3)
	out[2] = graphMatching->serialize_D2();
    delete graphMatching;
}
