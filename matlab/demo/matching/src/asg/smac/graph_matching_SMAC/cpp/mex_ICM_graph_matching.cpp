//function [X,scores]=mex_ICM_graph_matching(X,E12,Ws,Wd);
// Timothee Cour, 4 february 2009: fixed a bug which affected discretization in the case n1~=n2

//TODO:handle sparse Ws
#undef NDEBUG
#include <math.h>
#include <mex.h>
# include "mex_util.cpp"
#include "Matrix.cpp"
# include "mex_math.cpp"
#include <vector>

#include "mex_time.cpp"


double compute_delta_swap_quadratic(double*Ws,double*Wd,double*X,int*indexes,int i,int j,int n){
    //Ws:nxn,Wd:nx1,X:nx1
    double score=0.5*(Wd[j]-Wd[i]);
    double*pWs1=&Ws[i*n];
    double*pWs2=&Ws[j*n];
    for(int u=0;u<n;u++){
        if(u!=i && X[indexes[u]]){
            score+=*pWs2-*pWs1;
        }
        pWs1++;
        pWs2++;
    }
    return 2*score;    
}
double compute_delta_swap_quadratic_sparse(double*pr_Ws,mwIndex*ir_Ws,mwIndex*jc_Ws,double*Wd,double*X,int*indexes,int i,int j,int n){
    //Ws:nxn,Wd:nx1,X:nx1
    double score=0.5*(Wd[j]-Wd[i]);
    int v;
    for(v=jc_Ws[j];v<jc_Ws[j+1];v++){
        int u=ir_Ws[v];
        if(u!=i && X[indexes[u]]){
            score+=pr_Ws[v];
        }
    }
    for(v=jc_Ws[i];v<jc_Ws[i+1];v++){
        int u=ir_Ws[v];
        if(u!=i && X[indexes[u]]){
            score-=pr_Ws[v];
        }
    }
    return 2*score;    
}

class Class1{
    public:
        const mxArray *mxWs;
        Matrix<double>X;
        Matrix<char>E12;
//         Matrix<double>Ws;
        Matrix<double>Wd;
        Matrix<double>scores;
        Matrix<int>match_index;
        vector<int>indexes;
        Matrix<int>labels;
        vector<int>labels_free;
        int n,k;
        int nbMatches;
        double valmin;
//         double score0;
        void mexGate(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
            X.readmxArray(in[0]);
            E12.readmxArray(in[1]);
//             Ws.readmxArray(in[2]);
            mxWs=in[2];
            
            Wd.readmxArray(in[3]);
//             score0=mxGetScalar(in[4]);
            n=E12.p;
            k=E12.q;
            assert(n<=k);
            
            compute_match_index();
            compute_labels();
//             valmin=mxGetEps();
            valmin=1e-10;
//             valmin=1e-5;
            ICM_graph_matching();
            out[0]=X.tomxArray();
            out[1]=scores.tomxArray();
            
        }
        void compute_labels(){
            labels.resize(n,1);
            int i,a;
            for(i=0;i<n;i++){
                for(a=0;a<k;a++){
                    if(X[i+a*n]){
                        labels[i]=a;
                        break;
                    }
                }
            }
            
            labels_free.resize(0);
            for(a=0;a<k;a++){
                bool is_free=true;
                for(i=0;i<n;i++){
                    if(X[i+a*n]){
                        is_free=false;
                        break;
                    }
                }
                if(is_free)
                    labels_free.push_back(a);
            }
            assert(labels_free.size()==k-n);
        }
        void compute_match_index(){
            nbMatches=0;
            match_index.resize(n,k);
            for(int i=0;i<n*k;i++){
                if(E12[i]){
                    match_index[i]=nbMatches++;
                    indexes.push_back(i);
                }
            }
        }
        
        void ICM_graph_matching(){
            bool isChanged=true;
            int nbIter=0;
            while(isChanged){
                isChanged=false;
                double score=0;
                for(int i=0;i<n;i++){
                    ICM_graph_matching_node(i,isChanged,score);                    
                }
                scores.push_back(score);
                nbIter++;
                if(nbIter==1000){
                    Matrix<double>::disp(scores.data);
                }
                assert(nbIter<1000);//VOIR
            }
        }
        double compute_score_swap_match(int i,int a,int b){
            int u=match_index[i+a*n];
            int v=match_index[i+b*n];
            if(mxIsSparse(mxWs))
                return compute_delta_swap_quadratic_sparse(mxGetPr(mxWs),mxGetIr(mxWs),mxGetJc(mxWs),&Wd[0],&X[0],&indexes[0],u,v,nbMatches);
            else
                return compute_delta_swap_quadratic(mxGetPr(mxWs),&Wd[0],&X[0],&indexes[0],u,v,nbMatches);
        }
        double compute_score_swap_match_pair(int i,int a,int j,int b){
            double score1=compute_score_swap_match(i,a,b);
            X[i+a*n]=0;
            X[i+b*n]=1;
            double score2=compute_score_swap_match(j,b,a);
            X[i+b*n]=0;
            X[i+a*n]=1;
            return score1+score2;
        }

        void ICM_graph_matching_node(int i,bool&isChanged,double&score){
            int a=labels[i];
            assert(X[i+a*n]==1);           
            
            bool isValid=false;
            double score_best=0;
            int j_best=0;
            double score_i=0;
            for(int j=0;j<n;j++){
                int b=labels[j];
                assert(X[j+b*n]==1);
                
                if(!(E12[i+n*b] && E12[j+n*a]))
                    continue;
                double score_j=compute_score_swap_match_pair(i,a,j,b);
                if(j==i)
                    score_i=score_j;
                if(!isValid||score_j>score_best){
                    score_best=score_j;
                    isValid=true;
                    j_best=j;                    
                }                                
            }
            assert(isValid);
//             assert(fabs(score_i)<1e-10);
            
//             assert(score_best>=0);
//             if(!(score_best>=-10*valmin))
//                 DEBUG(score_best);
//             assert(score_best>=-10*valmin);
//             DEBUG(score_best);
            
            if(!labels_free.empty()){
                double score_change;
                int u_change;
                bool isValid_change;
                compute_score_best_change_free(i, a, score_change, u_change,isValid_change);
                if(isValid_change && (score_change > score_best && score_change > score_i+valmin)){
                    assert(score_change>0);
                    isChanged=true;
                    int b=labels_free[u_change];
                    labels_free[u_change]=a;
                    labels[i]=b;
                    X[i+a*n]=0;
                    X[i+b*n]=1;
                    score+=score_change;
                    return;
                }
            }
            
            if(j_best==i || score_best<=score_i+valmin){
                return;
            }
            
            assert(score_best>0);
            
            isChanged=true;
            int b=labels[j_best];
            labels[i]=b;
            labels[j_best]=a;
            X[i+a*n]=0;
            X[i+b*n]=1;
            X[j_best+b*n]=0;
            X[j_best+a*n]=1;                       
            score+=score_best;
        }        
        void compute_score_best_change_free(int i,int a,double&score,int&u_change,bool&isValid){
            isValid=false;
            for(int u=0;u<labels_free.size();u++){
                int bu=labels_free[u];
                if(!(E12[i+n*bu]))
                    continue;
                double scoreu=compute_score_swap_match(i,a,bu);
                if(!isValid || scoreu>score){
                    isValid=true;
                    score=scoreu;
                    u_change=u;
                }
            }
        }
};


void mexFunction(int nargout, mxArray *out[], int nargin, const mxArray *in[]) {
    Class1 class1;
    class1.mexGate(nargout,out,nargin,in);
}

