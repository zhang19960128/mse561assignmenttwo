#include "atom.h"
#include "parameter.h"
#include <cmath>
#include <vector>
#include "ndarrays.h"
double distance(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	return sqrt(r);
}
double potential(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	r=sqrt(r);
	if(r<r0){
		return 4*eps*(pow(sigma/r,12)-pow(sigma/r,6));
	}
	else if(r>r_cut){
		return 0;
	}
	else{
		return A*pow(r-r_cut,3)+B*pow(r-r_cut,2);
	}
}
//the force exerted on one by two
std::vector<double> str_tensor(atom& one,atom& two){
	std::vector<double> a(4,0);
	double r=distance(one,two);
    double deri=0;
    if(r<r0){
        deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
    }
    else if(r<r_cut){
        deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
    }
    else{
        deri=0.0;
    }
    a[0]=deri*(one.x-two.x)*(one.x-two.x)/r;
    a[1]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[2]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[3]=deri*(one.y-two.y)*(one.y-two.y)/r;
    return a;
}
void updatelist(ndarrays<atom>& input,int size){
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++){
			input(i,j).neighborx.clear();
			input(i,j).neighbory.clear();
		}
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++)
			for(size_t k=0;k<size;k++)
				for(size_t l=0;l<size;l++)
				{
					if(distance(input(i,j),input(k,l))<r_cut){
						input(i,j).neighborx.push_back(k);
						input(i,j).neighbory.push_back(l);
					}
				}
}