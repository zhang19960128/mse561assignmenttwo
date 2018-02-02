#ifndef atom_h
#define atom_h
#include "parameter.h"
#include <list>
#include <vector>
#include "ndarrays.h"
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
