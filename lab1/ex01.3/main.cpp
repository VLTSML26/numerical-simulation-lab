/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
 
int main() {

	Random rnd;

  	const int N_throws = 10000;
  	const int N_blocks = 100;
  	const double d = 1., l = 0.7;

  	ofstream write("ex01.3.out");
	for(int i=0; i<N_blocks; i++){
     	int hit = 0;
     	for(int j=0; j<N_throws; j++){
       		double x1 = rnd.Rannyu(0., d);
       		double x, y, theta;
       		do {
       			x = rnd.Rannyu();
       			y = rnd.Rannyu();
       		} while (x*x + y*y > 1.);
      		theta = acos(x/sqrt(x*x + y*y)); // NB: y sempre >= 0
       		double x2 = x1 + l*sin(theta);
       		if(x2 < 0. or x2 > d) hit++;
     	}
		double pi = 2*l*N_throws / (hit*d);
     	write << pi << endl;
	}	

	rnd.SaveSeed();
   	return 0;
}
