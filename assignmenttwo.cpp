#include<iostream>
#include "parameter.h"
#include "atom.h"
#include "ndarrays.h"
#include <list>
#include <fstream>
int main(){
	int size=20;
	ndarrays<atom> atomall(2,size,size);
	std::fstream atomfile;
	atomfile.open("atominfo.txt",std::fstream::out);
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++){
			atomall(i,j).setx(i*r_min);
			atomall(i,j).sety(j*r_min);
			atomall(i,j).setr(0.1*r_min);
			atomfile<<i*r_min<<" "<<j*r_min<<" "<<0.1*r_min<<std::endl;
		}
    atomfile.close();
	updatelist(atomall,size);
    updatetensor(atomall,size);
    double r_start=0.0000001;
    double r_stop=15;
    int N=10000;
    std::vector<double> ra_dis(N,0.0);
    double r_delta=(r_stop-r_start)/N;
    double r_inter=0.0;
    int count_old=0;
    int count_new=0;
    int count_delta=0;
    std::fstream radis;
    radis.open("radius_distribution.txt",std::fstream::out);
    for(size_t i=0;i<N;i++){
        r_inter=i*r_delta+r_start;
        count_new=count(atomall,atomall(10,10),r_inter,size);
        count_delta=count_new-count_old;
        ra_dis[i]=count_delta/2/Pi/r_inter/r_delta;
        count_old=count_new;
        radis<<r_inter<<" "<<ra_dis[i]<<std::endl;
    }
    radis.close();
}
