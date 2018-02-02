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
	updatelist(atomall,size);

}
