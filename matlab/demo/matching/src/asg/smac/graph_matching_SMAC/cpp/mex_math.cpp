/*================================================================
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

mex_math = used by a couple of mex functions 
 *=================================================================*/
#pragma once
# include "math.h"
#include <vector>
#include "Matrix.cpp"

int round2(double x) {
    //return floor(x+0.5);
    return x>=0 ? (int)(x+0.5) : (int) (x-0.5);
}
/* Problem: when compiling on opteron, says error: new declaration
int round(double x) {
    //return floor(x+0.5);
    return x>=0 ? (int)(x+0.5) : (int) (x-0.5);
}*/

double min(double x,double y) {
    return x<y?x:y;
}
double max(double x,double y) {
    return x>y?x:y;
}
int max(int x,int y) {
    return x>y?x:y;
}
int min(int x,int y) {
    return x<y?x:y;
}

double vec_max(double *x,int n) {
    double res=x[0];
    for(int i=1;i<n;i++)
        if(res<x[i])
            res=x[i];
    return res;
}
int vec_max(int *x,int n) {
    int res=x[0];
    for(int i=1;i<n;i++)
        if(res<x[i])
            res=x[i];
    return res;
}
int vec_min(int *x,int n) {
    int res=x[0];
    for(int i=1;i<n;i++)
        if(res>x[i])
            res=x[i];
    return res;
}
double vec_min(double *x,int n) {
    double res=x[0];
    for(int i=1;i<n;i++)
        if(res>x[i])
            res=x[i];
    return res;
}
double dot(double *x,double *y, int n) {
    double temp = 0;
    for(int i=0;i!=n;i++)
        temp+=*x++ * *y++;
    return temp;
}

void scalar_times_vec(double a, double *x,double *y, int n) {
    // y=a*x
    for(int i=0;i!=n;i++)
        *y++ = *x++ * a;
}
void scalar_plus_vec_self(int a, int *x, int n) {
    // x=a+x
    for(int i=0;i!=n;i++)
        *x++ += a;
}

void vec_plus_vec(double *x1, double *x2,double *y, int n) {
    // y=x1+x2
    for(int i=0;i!=n;i++)
        y[i] = x1[i] + x2[i];
        //*y++ = *x1++ + *x2++;
}


void ind2sub(int*ind,int*indi,int*indj,int p,int q,int n){
	//caution: output is in matlab conventions//TODO;correct this [mex_ind2sub
    int i;
    int*pindi=indi;
    int*pindj=indj;
    int*pind=ind;
    int temp;
    for(i=0;i<n;i++){
        temp=*pind++-1;
        *pindi++=temp%p+1;
        *pindj++=temp/p+1;
    }
}

//void symmetrize_sparse()

void set_enum(vector<int>&indexes,int x1,int x2,bool flip){
    int nb=x2-x1+1;
    assert(nb>=0);
    indexes.resize(nb);
    int temp=0;
    if(flip){
        temp=x2;
        for(int i=0;i<nb;i++){
            indexes[i]=temp--;
        }
    }
    else{
        temp=x1;
        for(int i=0;i<nb;i++){
            indexes[i]=temp++;
        }
    }
}

void intline2(int x1,int y1,int x2,int y2,vector<int>&xs,vector<int>&ys){
    int dx=abs(x2-x1);
    int dy=abs(y2-y1);
    if(dx==0 && dy==0){
        xs.push_back(x1);
        ys.push_back(y1);
        return;
    }
    bool flip=false;
    int t;
    if(dx>=dy){
        if(x1>x2){
            t=x1;x1=x2;x2=t;
            t=y1;y1=y2;y2=t;
            flip=true;
        }
        double m=((double)(y2-y1)) / ((double)(x2-x1));
        set_enum(xs,x1,x2,flip);
        ys.resize(xs.size());
        for(int u=0;u<xs.size();u++){
            ys[u]=0.5+m*(xs[u]-x1)+y1;
        }
//         Matrix<int>::disp(xs);
    }
    else{
        if(y1>y2){
            t = x1; x1 = x2; x2 = t;
            t = y1; y1 = y2; y2 = t;
            flip = true;
        }
        double m=((double)(x2-x1)) / ((double)(y2-y1));
        set_enum(ys,y1,y2,flip);
        xs.resize(ys.size());
        for(int u=0;u<ys.size();u++){
            xs[u]=0.5+m*(ys[u]-y1)+x1;
        }
//         Matrix<int>::disp(ys);
    }
}
