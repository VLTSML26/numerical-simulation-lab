/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_throws	10000
#define N_blocks	100
#define d			1.
#define l			0.7

#include "../../Random/Random.h"
#include <fstream>
#include <cmath>

using namespace std;
 
int main() {

	// initialize random and output
	Random rnd;
  	ofstream write("data/ex01.3.out");

	// data blocking method for pi with Buffon's experiment
	for(int i=0; i<N_blocks; ++i){
     	int hit = 0;
     	for(int j=0; j<N_throws; ++j){

			// sets first endpoint of the needle
       		double x1 = rnd.Rannyu(0., d);

			// generates a random angle theta WITHOUT USING pi
			// (note that y is always positive)
       		double x, y, theta;
       		do {
       			x = rnd.Rannyu();
       			y = rnd.Rannyu();
       		} while (x*x + y*y > 1.);
      		theta = acos(x/sqrt(x*x + y*y));

			// sets second endpoint of the needle
       		double x2 = x1 + l*sin(theta);

			// checks if it hit
       		if(x2 < 0. or x2 > d) hit++;
     	}
		
		// find pi from Buffon's formula
		double pi = 2*l*N_throws / (hit*d);
     	write << pi << endl;
	}	

	rnd.SaveSeed();
   	return 0;
}
