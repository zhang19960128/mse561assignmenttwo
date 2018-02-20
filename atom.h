#ifndef atom_h
#define atom_h
#include "parameter.h"
#include <list>
#include <vector>
#include "ndarrays.h"
#include <iostream>
#include <fstream>
class atom{
	public:
		atom()=default;
		atom(double x1,double x2,double r):x(x1),y(x2),radius(r){};
		friend double distance(atom&,atom&);
		friend double potential(atom&,atom&);
		friend void updatelist(ndarrays<atom>&,int);
    friend void updatetensor(ndarrays<atom>&,int);
		friend std::vector<double> str_tensor(atom&,atom&);
    void setx(double x1){
			x=x1;
		}
		void sety(double x2){
			y=x2;
		}
		void setr(double ra){
			radius=ra;
		}
    void setstress_tensor(std::vector<double> a){
            stresstensor=a;
    }
    std::vector<double> getstress(){
			return stresstensor;
    }
		double getx(){
			return x;
		}
		double gety(){
			return y;
		}
		double getr(){
			return radius;
		}
	void printneighbor(){
			std::list<int>::iterator b=neighbory.begin();
			std::cout<<"the neighbors are:"<<std::endl;
			for(std::list<int>::iterator a=neighborx.begin();a!=neighborx.end();a++){
				std::cout<<"("<<*a<<","<<*b<<")"<<std::endl;
				b++;
			}
		}
	void printstress(){
		std::cout<<"the stress tensor is: (sigma_xx,sigma_xy;sigma_yx,sigma_yy):"<<std::endl;
		std::cout<<stresstensor[0]<<" "<<stresstensor[1]<<std::endl;
		std::cout<<stresstensor[2]<<" "<<stresstensor[3]<<std::endl;
	}
	void printinfo(){
		std::cout<<"for atom position ("<<x<<","<<y<<")"<<std::endl;
		printneighbor();
		printstress();
	}
	friend std::ostream& operator<<(std::ostream& os,atom& output);
	friend std::fstream& operator<<(std::fstream& fs,atom& output);
	private:
		double x;
		double y;
		double radius;
		std::list<int> neighborx;
		std::list<int> neighbory;
    std::vector<double> stresstensor;
};
int count(ndarrays<atom>& all,atom&,double r,int size);
#endif
