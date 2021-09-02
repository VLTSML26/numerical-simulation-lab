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

#include "../../Random/Random.h"
#include <fstream>
#include <cmath>

using namespace std;
 
int main() {

	// initialize random and outputs
	Random rnd;
  	ofstream unif("data/ex02.1_uniforme.out");
  	ofstream linear("data/ex02.1_lineare.out");

	// data blocking method with uniform sampling
	for(int i=0; i<N_blocks; i++) {
     	double sum = 0.;
     	for(int j=0; j<N_throws; j++) {
       		double ran = rnd.Rannyu();
       		double y = M_PI / 2 * cos(M_PI / 2 * ran);
			sum += y;
     	}
		sum /= N_throws;
     	unif << sum << endl;
	}
	
	// data blocking method with importance sampling
	for(int i=0; i<N_blocks; i++) {
     	double sum = 0.;
     	for(int j=0; j<N_throws; j++) { 

			// p(x) = 2(1-x)
       		double ran = rnd.Line();
       		double y = (M_PI / 2 * cos(M_PI / 2 * ran)) / (2 * (1 - ran));
			sum += y;
     	}
		sum /= N_throws;
     	linear << sum << endl;
	}
	
	rnd.SaveSeed();
   	return 0;
}
