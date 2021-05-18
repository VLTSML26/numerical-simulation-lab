/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "../../Random/Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
 
int main() {

	Random rnd;

  	const int N_throws = 10000;
  	const int N_blocks = 100;
	
  	ofstream write("ex02.1.1.out"); // uniforme
	for(int i=0; i<N_blocks; i++) {
     	double sum = 0.;
     	for(int j=0; j<N_throws; j++) {
       		double ran = rnd.Rannyu();
       		double y = M_PI / 2 * cos(M_PI / 2 * ran);
			sum += y;
     	}
		sum /= N_throws;
     	write << sum << endl;
	}
	write.close();
	
  	ofstream write2("ex02.1.2.out"); // p(x) = 2(1-x)
	for(int i=0; i<N_blocks; i++) {
     	double sum = 0.;
     	for(int j=0; j<N_throws; j++) {
       		double ran = rnd.Line();
       		double y = (M_PI / 2 * cos(M_PI / 2 * ran)) / (2 * (1 - ran));
			sum += y;
     	}
		sum /= N_throws;
     	write2 << sum << endl;
	}
	write2.close();
	
	rnd.SaveSeed();
   	return 0;
}
  	/*ofstream write1("ex02.1.out"); 
	for(int i=0; i<N_blocks; i++) {
     	double sum[2] = {0.};
     	for(int j=0; j<N_throws; j++) {
       		double ran[2] = {rnd.Rannyu(), rnd.Line()};
       		double y[2] = {0.};
			y[0] = M_PI / 2 * cos(M_PI / 2 * ran[0]);
			y[1] = (M_PI / 2 * cos(M_PI / 2 * ran[1])) / (2 * (1 - ran[1]));
			for(int k=0; k<2; k++)
				sum[k] += y[k];
     	}
		for(int k=0; k<2; k++)
			sum[k] /= N_throws;
     	write1 << sum[0] << "\t" << sum[1] <<  endl;
	}
	write1.close();*/
