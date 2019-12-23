/*
// Timothee Cour, 21-Apr-2008 17:31:23
// This software is made publicly for research use only.
// It may be modified and redistributed under the terms of the GNU General Public License.

TODO: do a macro to enable/disable timing

usage:
#include "mex_time.cpp"
resetTimer();
dispTimer();
tic(1);
toc(1);

*/

#pragma once
#include <time.h>
#include <map>
#include <sstream>

//void tic(int s);
//void toc(int s);
void tic(string s);
void toc(string s);
string int2string(int i);

//# define TICTOC( expr ) tic((int)__LINE__); expr ;toc((int)__LINE__)
//# define TICTOC( expr ) tic(string(#expr)); expr ;toc(string(#expr))
//# define TICTOC( expr ) tic(string(#expr)+(char*)__LINE__); expr ;toc(string(#expr)+(char*)__LINE__)
//# define TICTOC( expr ) tic(string(#expr)+"\""__LINE__"\""); expr ;toc(string(#expr)+"\""__LINE__"\"")
# define TICTOC( expr ) tic(int2string(__LINE__)+", "+string(#expr)); expr ;toc(int2string(__LINE__)+", "+string(#expr))


string int2string(int i){
	string s;
	ostringstream oss;
	oss << i;
	return oss.str();
}


long timeCounter;
map<string,long>timeCounter_map;
map<string,int>timeCounter_map_nb;
map<string,long>timeCounter_map_aux;

vector<long>timeCounter_vec;
vector<int>timeCounter_vec_nb;
vector<bool>timeCounter_vec_parity;
vector<long>timeCounter_vec_aux;

bool is_dispTimer=true;

void tic(){
	timeCounter=clock();
}

void toc(){
	long elapsed_time=clock()-timeCounter;
	mexPrintf("Elapsed time (in C) is = %10.5g\n",((double)elapsed_time)/CLOCKS_PER_SEC);   
}
void tic(string s){
	timeCounter_map_aux[s]=clock();
}
void toc(string s){
	timeCounter_map[s]+=clock()-timeCounter_map_aux[s];
	timeCounter_map_nb[s]++;
}
void tic(int s){
	if(s>=timeCounter_vec.size()){
		timeCounter_vec.resize(s+1);
		timeCounter_vec_aux.resize(s+1);
		timeCounter_vec_nb.resize(s+1);
	}
	timeCounter_vec_aux[s]=clock();
}
void tic2(int s){
	if(s>=timeCounter_vec.size()){
		timeCounter_vec.resize(s+1);
		timeCounter_vec_aux.resize(s+1);
		timeCounter_vec_nb.resize(s+1);
		timeCounter_vec_parity.resize(s+1);
	}
	timeCounter_vec_parity[s]=!timeCounter_vec_parity[s];
	if(timeCounter_vec_parity[s])
		timeCounter_vec_aux[s]=clock();
	else{
		timeCounter_vec[s]+=clock()-timeCounter_vec_aux[s];
		timeCounter_vec_nb[s]++;
	}

}

void toc(int s){
	assert(s<timeCounter_vec.size());
	assert(s<timeCounter_vec_nb.size());
	assert(s<timeCounter_vec_aux.size());
	timeCounter_vec[s]+=clock()-timeCounter_vec_aux[s];
	timeCounter_vec_nb[s]++;
}

void resetTimer(){
	timeCounter_map.clear();
	timeCounter_map_aux.clear();
	timeCounter_map_nb.clear();

	timeCounter_vec.clear();
	timeCounter_vec_aux.clear();
	timeCounter_vec_nb.clear();
	timeCounter_vec_parity.clear();

	timeCounter=0;
}

void dispTimer(){

	if(!is_dispTimer)
		return;


	double sum=0,timei;

	int maxSize=0;
	int maxSizei;
	map<string,long>::iterator iter;
	for(iter = timeCounter_map.begin(); iter != timeCounter_map.end(); iter++){
		maxSizei=(*iter).first.size();
		if(maxSize<maxSizei)
			maxSize=maxSizei;
	}
	/*
	for(iter = timeCounter_map.begin(); iter != timeCounter_map.end(); iter++){
		timei=((double)(*iter).second)/CLOCKS_PER_SEC;
		mexPrintf("timing[%s] = %10.5g\n",(*iter).first.c_str(),timei);   
		sum+=timei;
    }
	*/

	map<string,int>::iterator iter2=timeCounter_map_nb.begin();
	for(iter = timeCounter_map.begin(); iter != timeCounter_map.end(); iter++,iter2++){
		timei=((double)(*iter).second)/CLOCKS_PER_SEC;
		//mexPrintf("timing[%s] = %10.5g\n",(*iter).first.c_str(),timei);   
		int nb=(*iter2).second;
		double timeiPerIter = timei/nb;
		mexPrintf("timing[\"%-*s\"] = %10.3g   %10.3g iter * %10.3g /iter\n",maxSize,(*iter).first.c_str(),timei,(double)nb,timeiPerIter);   
    }


	for(int i=0;i<timeCounter_vec.size();i++){
		int nb=timeCounter_vec_nb[i];
		if(nb==0)
			continue;
		timei=((double)timeCounter_vec[i])/CLOCKS_PER_SEC;
		double timeiPerIter = timei/nb;
		mexPrintf("timing[%-*d] = %10.3g   %10.3g iter * %10.3g /iter\n",maxSize+2,i,timei,(double)nb,timeiPerIter);   
	}
}