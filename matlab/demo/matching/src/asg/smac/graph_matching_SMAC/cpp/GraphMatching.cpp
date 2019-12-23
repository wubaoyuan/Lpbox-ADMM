// options.cliqueModeString = {"log_linear", "L1", "edges", "binary"}
// options.normalization = {"none", "D1", "D2", "D1D2", "iterative"}
// options.nbIter = int
// f.fG = double[size(G, 3)]
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <float.h>
#include <algorithm>
#include "mex_util.cpp"
#include "mex_math.cpp"

class GraphMatching {
    mxArray *W;
    const mxArray *mxArray_W1;
    const mxArray *mxArray_W2;
    double *G1;
    double *G2;
    double *W12;
    vector<int> I12;
    double *F12;
    double *fG;
    double *valW;

    int n1, n2, n12;
    int nnz;
    int RG, RF;

    vector<int> indj;
    mwIndex* indi;
    mwIndex* jc;

    vector<double> D1;
    vector<double> D2;

    string normalization;
    int nbIter;

    enum CliqueMode {
	mode_L1, mode_log_linear, mode_edges, mode_binary, mode_additive
    };
    string cliqueModeString;
    CliqueMode cliqueMode;

public: 

    GraphMatching() {}

    ~GraphMatching() {}

    void compute_W() {
	nnz = compute_nnz_W();
	W = mxCreateSparse(n12, n12, nnz, mxREAL);
	valW = mxGetPr(W);

	indj.resize(nnz);

	indi = mxGetIr(W); // voir: uniformiser (vector<int>)?
	compute_ind_W();

	jc=mxGetJc(W);
	compute_j2Jc(jc, &indj[0], n12, nnz);
	compute_valWs();
    }

    void normalize_W() {
	if (normalization == "none") {
	    ;
	} else if (normalization == "D1") {
	    computeD1();
	    normalize_D1();
	} else if (normalization == "D2") {
	    computeD2();
	    normalize_D2();
	} else if (normalization == "D1D2"){
	    computeD1();
	    computeD2();
	    normalize_D1D2();
	} else if (normalization == "iterative") {
	    normalize_iterative(nbIter);
	} else {
	    printError;
	}
    }

    void deserializeGraphs(const mxArray *mxArray_W1,
			   const mxArray *mxArray_W2,
			   const mxArray *mxArray_I12,
			   const mxArray *mxArray_W12,
			   const mxArray *mxArray_G1,
			   const mxArray *mxArray_G2,
			   const mxArray *mxArray_F12,
			   const mxArray *mxArray_f) {
	this->mxArray_W1 = mxArray_W1;
	this->mxArray_W2 = mxArray_W2;
	this->G1 = mxGetPr(mxArray_G1);
	this->G2 = mxGetPr(mxArray_G2);			
	this->W12 = mxGetPr(mxArray_W12);
	this->F12 = mxGetPr(mxArray_F12);
	this->fG = mxGetPr(getFieldByName(mxArray_f, 0, "fG"));
	
	// this->fF = mxGetPr(getFieldByName(mxArray_f, 0, "fF"));
	n1 = mxGetM(mxArray_W1);
	n2 = mxGetM(mxArray_W2);

	// assert(n1 == mxGetM(mxArray_W1) && n2 == mxGetM(mxArray_W2));
	this->RG = getSize3(mxArray_G1);
	this->RF = getSize3(mxArray_F12);
	deserialize_I12(mxArray_I12);
    }

    void deserializeOptions(const mxArray* options) {
	if (isField(options, "n1n2")) {
	    vector<double> n1n2;
	    mxArrayStruct2array(options, 0, "n1n2", n1n2, 0);
	    assert(n1n2.size() == 2);
	    n1 = n1n2[0];
	    n2 = n1n2[1];
	}

	if (isField(options, "normalization"))
	    mxArrayStruct2array(options, 0, "normalization", normalization);

	if (normalization=="iterative")
	    assert(isField(options,"nbIter"));
	if (isField(options,"nbIter")) {
	    mxArrayStruct2array(options,0,"nbIter",nbIter,0);
	}
	if (isField(options,"cliqueModeString")) {
	    mxArrayStruct2array(options, 0,"cliqueModeString",cliqueModeString);
	    if(cliqueModeString=="log_linear")
		cliqueMode=mode_log_linear;
	    else if(cliqueModeString=="L1")
		cliqueMode=mode_L1;
	    else if(cliqueModeString=="edges")
		cliqueMode=mode_edges;
	    else if(cliqueModeString=="binary")
		cliqueMode=mode_binary;
	    else if(cliqueModeString=="additive")
		cliqueMode=mode_additive;
	    else
		printError;
	} else{
	    cliqueModeString="log_linear";
	    cliqueMode=mode_log_linear;
	}
    }

    void deserialize_I12(const mxArray *mxArray_I12){
	mxArray2array(mxArray_I12, I12, -1);		
	n12 = mxGetM(mxArray_I12);
	int n1_max = *max_element(I12.begin(), I12.begin() + n12) + 1;
	int n2_max = *max_element(I12.begin() + n12, I12.end()) + 1;
	assert(n1_max <= n1 && n2_max <= n2);
    }

    void deserializeW(const mxArray*W0, const mxArray *mxArray_I12) {
	deserialize_I12(mxArray_I12);
	assert(mxGetM(W0) == n12);
	W = mxDuplicateArray(W0);
	valW = mxGetPr(W);
	indi = mxGetIr(W);
	jc = mxGetJc(W);
	nnz = jc[n12];
    }

    int compute_nnz_W() {
	mwIndex *ir1 = mxGetIr(mxArray_W1);
	mwIndex *jc1 = mxGetJc(mxArray_W1);
	int n1 = mxGetM(mxArray_W1);
	mwIndex *ir2 = mxGetIr(mxArray_W2);
	mwIndex *jc2 = mxGetJc(mxArray_W2);
	int n2 = mxGetM(mxArray_W2);

	int temp;
	int nnz = 0;
	int k, i1, j1, i2, j2, k2, k1;

	for (k = 0; k < n12; k++) {
	    j1 = I12[k];
	    j2 = I12[k + n12];
	    for (k2 = jc2[j2]; k2 != jc2[j2 + 1]; k2++) {
		i2 = ir2[k2];
		temp = n1 * i2;

		for (k1 = jc1[j1]; k1 != jc1[j1 + 1]; k1++) {
		    i1 = ir1[k1];
		    if (W12[i1 + temp] && (i2 > j2 || (i2 == j2 && i1 >= j1)))
			//to keep lower part of W
			nnz++;
		}
	    }
	}
	return nnz;
    }

    void compute_ind_W() {
	mwIndex *ir1 = mxGetIr(mxArray_W1);
	mwIndex *jc1 = mxGetJc(mxArray_W1);
	int n1 = mxGetM(mxArray_W1);
	mwIndex *ir2 = mxGetIr(mxArray_W2);
	mwIndex *jc2 = mxGetJc(mxArray_W2);
	int n2 = mxGetM(mxArray_W2);
	int temp;

	int k,i1,j1,i2,j2,k2,k1;
	int nnz = 0;
	for(k=0;k<n12;k++) {
	    j1 = I12[k];
	    j2 = I12[k+n12];
	    for (k2=jc2[j2]; k2!=jc2[j2+1]; k2++) {
		i2 = ir2[k2];
		temp = n1*i2;
		for (k1=jc1[j1]; k1!=jc1[j1+1]; k1++) {
		    i1 = ir1[k1];
		    if (W12[i1+temp] && (i2>j2 || (i2==j2&&i1>=j1) )) {//to keep lower part of W
			indi[nnz] = W12[i1+temp]-1;
			indj[nnz] = k;
			nnz++;
		    }
		}
	    }
	}
    }

    void compute_valWs() {
	int i1,i2,j1,j2,i;
	for(int j=0;j<n12;j++){
	    for(int k=jc[j];k<jc[j+1];k++){
		i=indi[k];
		i1=I12[i];
		i2=I12[i+n12];
		j1=I12[j];				
		j2=I12[j+n12];
		valW[k] = computeClique(i1,i2,j1,j2);
	    }
	}
    }

    void computeD1() {
	//D1.resize(n1*n1,DBL_EPSILON);
	D1.resize(n1 * n1);
	D1.assign(n1 * n1, DBL_EPSILON);
	double val;
	int i1, i2, j1, j2, i;
	for (int j = 0; j < n12; j++) {
	    for (int k = jc[j]; k < jc[j + 1]; k++) {
		i = indi[k];
		i1 = I12[i];
		i2 = I12[i + n12];
		j1 = I12[j];				
		j2 = I12[j + n12];
		val = valW[k];
		D1[i1 + n1 * j1] += val >= 0 ? val : -val;
		if(i1 == j1 && i2 == j2)
		    continue; //don't count diagonal elements twice
		D1[j1 + n1 * i1] += val >= 0 ? val : -val; //because tril(W)
	    }
	}
    }

    void computeD2() {
	//D2.resize(n2*n2,DBL_EPSILON);
	D2.resize(n2 * n2);
	D2.assign(n2 * n2, DBL_EPSILON );
	double val;
	int i1, i2, j1, j2, i;
	for (int j = 0; j < n12; j++) {
	    for (int k = jc[j]; k < jc[j + 1]; k++) {
		i = indi[k];
		i1 = I12[i];
		i2 = I12[i + n12];
		j1 = I12[j];				
		j2 = I12[j + n12];
		val = valW[k];
		D2[i2 + n2 * j2] += val >= 0 ? val : -val;
		if(i1 == j1 && i2 == j2)
		    continue;//don't count diagonal elements twice
		D2[j2 + n2 * i2] += val >= 0 ? val : -val;//because tril(W)
	    }
	}
    }

    void normalize_D1() {
	int i1, i2, j1, j2, i;
	for (int j = 0; j < n12; j++) {
	    for (int k = jc[j]; k < jc[j + 1]; k++) {
		i = indi[k];
		i1 = I12[i];
		i2 = I12[i + n12];
		j1 = I12[j];				
		j2 = I12[j + n12];
		valW[k] /= D1[i1 + n1 * j1];				
	    }
	}
    }

    void normalize_D2() {
	int i1, i2, j1, j2, i;
	for (int j = 0; j < n12; j++){
	    for (int k = jc[j]; k < jc[j + 1]; k++) {
		i = indi[k];
		i1 = I12[i];
		i2 = I12[i + n12];
		j1 = I12[j];				
		j2 = I12[j + n12];
		valW[k] /= D2[i2 + n2 * j2];				
	    }
	}
    }

    void normalize_D1D2() {
	int i1,i2,j1,j2,i;
	for(int j=0;j<n12;j++){
	    for(int k=jc[j];k<jc[j+1];k++){
		i=indi[k];
		i1=I12[i];
		i2=I12[i+n12];
		j1=I12[j];				
		j2=I12[j+n12];
		valW[k]/=(D1[i1+n1*j1]*D2[i2+n2*j2]);		
	    }
	}
    }

    void normalize_iterative(int nbIter) {
	for(int i = 0; i < nbIter; i++) {
	    computeD1();
	    normalize_D1();
	    computeD2();
	    normalize_D2();
	}
    }

    virtual double computeClique(int i1,int i2,int j1,int j2){
	if(cliqueMode==mode_log_linear)
	    return computeClique_log_linear(i1,i2,j1,j2);
	if(cliqueMode==mode_L1)
	    return computeClique_L1(i1,i2,j1,j2);
	if(cliqueMode==mode_edges)
	    return computeClique_edges(i1,i2,j1,j2);
	if(cliqueMode==mode_binary)
	    return computeClique_binary(i1,i2,j1,j2);
	if(cliqueMode==mode_additive)
	    return computeClique_additive(i1,i2,j1,j2);
	printError;
    }
    
    virtual double computeClique_log_linear(int i1,int i2,int j1,int j2){
	double val=0;
	double w1,w2,w1T,w2T;
	int rG,rF;
	double m_i=1,m_j=1;//matching score correspondances for (i1,i2) and (j1,j2)

	for (rF=0;rF<RF;rF++) {
	    m_i*=F12[i1+n1*i2+n1*n2*rF];//voir
	    m_j*=F12[j1+n1*j2+n1*n2*rF];//voir
	}

	if(i1==j1 && i2==j2) {
	    //node-node similarity between i1 and i2
	    val = m_i*m_j;
	}
	else if(i1==j1 && i2!=j2 || i1!=j1 && i2==j2) {
	    //2 nodes are matched to 1 node
	    val = 0;
	    //val = 0.1;
	    // TODO : val = (score_of_merging i1 and j1 [resp,i2,j2] ) * m_i*m_j;
	    // val = m_i*m_j;
	}   
	else {
	    //make sure W(e1,e2)=W(e1T,e2T) ;
	    //but does not require W(e1,e2) = W(e2,e1)
	    for (rG=0;rG<RG;rG++) {
		w1 = G1[i1+n1*j1+n1*n1*rG];
		w2 = G2[i2+n2*j2+n2*n2*rG];
		w1T = G1[j1+n1*i1+n1*n1*rG];
		w2T = G2[j2+n2*i2+n2*n2*rG];
		//val+=sqrt(w1*w2)+sqrt(w1T*w2T);
		//val+=fG[rG]*(fabs(w1-w2)+fabs(w1T-w2T));//voir ;
		val+=fG[rG]*((w1-w2)*(w1-w2)+(w1T-w2T)*(w1T-w2T));//voir ;
	    }

	    //mexPrintf("i1,i2 j1,j2  w1,w2 = %d,%d %d,%d  %1.4g,%1.4g\n",i1,i2,j1,j2,w1,w2);
	    //val = exp(-val)*m_i*m_j;
	    val = exp(-val)*m_i*m_j;
	    //val = (exp(-val)-0.1)*m_i*m_j;
	}
	return val;
    }
    
    virtual double computeClique_additive(int i1,int i2,int j1,int j2){
	double val=0;
	double w1,w2,w1T,w2T;
	int rG,rF;
	double m_i=1,m_j=1;//matching score correspondances for (i1,i2) and (j1,j2)

	for (rF=0;rF<RF;rF++) {
	    m_i*=F12[i1+n1*i2+n1*n2*rF];//voir
	    m_j*=F12[j1+n1*j2+n1*n2*rF];//voir
	}

	if(i1==j1 && i2==j2) {
	    val = m_i*m_j;
	}
	else if(i1==j1 && i2!=j2 || i1!=j1 && i2==j2) {
	    val = 0;
	}   
	else {
	    for (rG=0;rG<RG;rG++) {
		w1 = G1[i1+n1*j1+n1*n1*rG];
		w2 = G2[i2+n2*j2+n2*n2*rG];
		w1T = G1[j1+n1*i1+n1*n1*rG];
		w2T = G2[j2+n2*i2+n2*n2*rG];
		if((w1==w2) || (w1T==w2T))
		    //if(w1==w2)
		    val+=fG[rG];
	    }
	    val = val*m_i*m_j;
	}
	return val;
    }

    virtual double computeClique_L1(int i1,int i2,int j1,int j2) {
	double val,w1,w2,w1T,w2T;
	int rG,rF;

	double m_i,m_j;//matching score correspondances for (i1,i2) and (j1,j2)
	//double temp;

	val = 0;

	m_i = 1;
	m_j = 1;
	for (rF=0;rF<RF;rF++) {
	    m_i*=F12[i1+n1*i2+n1*n2*rF];//voir
	    m_j*=F12[j1+n1*j2+n1*n2*rF];//voir
	}

	if(i1==j1 && i2==j2) {
	    //node-node similarity between i1 and i2
	    val = m_i*m_j;
	}
	else if(i1==j1 && i2!=j2 || i1!=j1 && i2==j2) {
	    //2 nodes are matched to 1 node
	    val = 0;
	    //val = 0.1;
	    // TODO : val = (score_of_merging i1 and j1 [resp,i2,j2] ) * m_i*m_j;
	    // val = m_i*m_j;
	}   
	else {
	    //make sure W(e1,e2)=W(e1T,e2T) ;
	    //but does not require W(e1,e2) = W(e2,e1)
	    for (rG=0;rG<RG;rG++) {
		w1 = G1[i1+n1*j1+n1*n1*rG];
		w2 = G2[i2+n2*j2+n2*n2*rG];
		w1T = G1[j1+n1*i1+n1*n1*rG];
		w2T = G2[j2+n2*i2+n2*n2*rG];
		val+=exp(-fG[rG]*(fabs(w1-w2)+fabs(w1T-w2T)));//voir ;
		//temp = 4.5-fabs(w1-w2)*fabs(w1-w2)/(2*0.5*0.5);
		//if(temp>0)
		//  val+=temp;
	    }
	    //mexPrintf("i1,i2 j1,j2  w1,w2 = %d,%d %d,%d  %1.4g,%1.4g\n",i1,i2,j1,j2,w1,w2);
	    val = val*m_i*m_j;
	}
	return val;
    }

    virtual double computeClique_edges(int i1,int i2,int j1,int j2) {
	double val,w1,w2,w1T,w2T;
	int rG,rF;

	double m_i,m_j;//matching score correspondances for (i1,i2) and (j1,j2)
	double temp;

	val = 0;

	m_i = 1;
	m_j = 1;
	for (rF=0;rF<RF;rF++) {
	    m_i*=F12[i1+n1*i2+n1*n2*rF];//voir
	    m_j*=F12[j1+n1*j2+n1*n2*rF];//voir
	}

	if(i1==j1 && i2==j2) {
	    //node-node similarity between i1 and i2
	    val = m_i*m_j;
	}
	else if(i1==j1 && i2!=j2 || i1!=j1 && i2==j2) {
	    //2 nodes are matched to 1 node
	    val = 0;
	    //val = 0.1;
	    // TODO : val = (score_of_merging i1 and j1 [resp,i2,j2] ) * m_i*m_j;
	    // val = m_i*m_j;
	}   
	else {
	    //make sure W(e1,e2)=W(e1T,e2T) ;
	    //but does not require W(e1,e2) = W(e2,e1)
	    for (rG=0;rG<RG;rG++) {
		w1 = G1[i1+n1*j1+n1*n1*rG];
		w2 = G2[i2+n2*j2+n2*n2*rG];
		w1T = G1[j1+n1*i1+n1*n1*rG];
		w2T = G2[j2+n2*i2+n2*n2*rG];

		temp = 0;
		if(rG==0)
		    {
			//temp = (w1+0.01)/(w2+0.01)-1;
			temp = 2*(w1-w2)/(w1+w2+0.0001);                    
			//val = 4.5-temp*temp/(2*0.2*0.2);
			temp = temp*temp;
		    }
		else if(rG==1)
		    {
			//temp = (1+cos(w1-w2))/2+(1+cos(w1T-w2T))/2;
			temp = (1-cos(w1-w2))*(1-cos(w1-w2));
		    }
		else if(rG==3)
		    //temp = (1+cos(w1-w2))/2+(1+cos(w1T-w2T))/2;
		    temp = (1-cos(w1-w2))*(1-cos(w1-w2));

		//temp = 4.5-fabs(w1-w2)*fabs(w1-w2)/(2*0.5*0.5);
		//val+=exp(-fG[rG]*(fabs(w1-w2)+fabs(w1T-w2T)));//voir ;
		//temp = 4.5-fabs(w1-w2)*fabs(w1-w2)/(2*0.5*0.5);
		val += fG[rG]*temp;
	    }
	    //mexPrintf("i1,i2 j1,j2  w1,w2 = %d,%d %d,%d  %1.4g,%1.4g\n",i1,i2,j1,j2,w1,w2);
	    //val = val*m_i*m_j;
	    val = exp(-val)*m_i*m_j;
	}
	return val;
    }


    virtual double computeClique_binary(int i1,int i2,int j1,int j2) {
	double val,w1,w2,w1T,w2T;
	int rG,rF;

	double m_i,m_j;//matching score correspondances for (i1,i2) and (j1,j2)
	val = 0;

	m_i = 1;
	m_j = 1;
	for (rF=0;rF<RF;rF++) {
	    m_i*=F12[i1+n1*i2+n1*n2*rF];//voir
	    m_j*=F12[j1+n1*j2+n1*n2*rF];//voir
	}

	if(i1==j1 && i2==j2) {
	    val = m_i*m_j;
	}
	else if(i1==j1 && i2!=j2 || i1!=j1 && i2==j2) {
	    //val = 0;
	    val = 0.1;
	}   
	else {
	    val = 1;
	    for (rG=0;rG<RG;rG++) {
		w1 = G1[i1+n1*j1+n1*n1*rG];
		w2 = G2[i2+n2*j2+n2*n2*rG];
		w1T = G1[j1+n1*i1+n1*n1*rG];
		w2T = G2[j2+n2*i2+n2*n2*rG];
		if( fabs(w1-w2)>=fG[rG] || fabs(w1T-w2T)>=fG[rG])
		    //val = 0;
		    val = -0.1;
	    }
	    val = val*m_i*m_j;
	}
	return val;
    }


	
    mxArray*serialize_indi1i2j1j2(){
	//n12,nnz,indi,jc
	mxArray* A=createMxArrayInt(nnz,4);
	int*pA=(int*)mxGetData(A);
	int i1,i2,j1,j2,i;
	for(int j=0;j<n12;j++){
	    for(int k=jc[j];k<jc[j+1];k++){
		i=indi[k];
		i1=I12[i];
		i2=I12[i+n12];
		j1=I12[j];				
		j2=I12[j+n12];
		pA[k+0*nnz]=i1+1;
		pA[k+1*nnz]=i2+1;
		pA[k+2*nnz]=j1+1;
		pA[k+3*nnz]=j2+1;
	    }
	}
	return A;
    }
    mxArray*serialize_indiW_indjW(){
	//n12,nnz,indi,jc
	mxArray* A=createMxArrayInt(nnz,2);
	int*pA=(int*)mxGetData(A);
	int i;
	for(int j=0;j<n12;j++){
	    for(int k=jc[j];k<jc[j+1];k++){
		i=indi[k];
		pA[k+0*nnz]=i+1;
		pA[k+1*nnz]=j+1;
	    }
	}
	return A;
    }

    mxArray* serialize_W() {
	return W;
    }
	
    mxArray* serialize_D1() {
	computeD1();
	//if(D1.size()==0)
	//	printError;
	mxArray *A = array2mxArray(D1);
	setDimensions(A, n1, n1);
	return A;
    }

    mxArray* serialize_D2() {	
	computeD2();
	//if(D2.size()==0)
	//	printError;
	mxArray* A = array2mxArray(D2);
	setDimensions(A, n2, n2);
	return A;
    }
};
