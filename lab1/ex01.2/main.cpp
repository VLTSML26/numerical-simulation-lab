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
  	const int Ns[4] = {1, 2, 10, 100};
  	const char stampa[4] = {'\t', '\t', '\t', '\n'};

  	ofstream write1 ("ex01.2.1.out"); // uniform
  	ofstream write2 ("ex01.2.2.out"); // exponential
  	ofstream write3 ("ex01.2.3.out"); // lorentz

	// ==== NOTA BENE ====
	// U sta per Uniform
	// E sta per Exponential
	// L sta per Lorentz
	
   	for(int i=0; i<N_throws; i++) {
     	for(int j=0; j<4; j++) {
       		double U_sum = 0.;
       		double E_sum = 0.;
       		double L_sum = 0.;
       		for(int k=0; k<Ns[j]; k++) {
       			U_sum += (int) rnd.Rannyu(1., 7.);
       			E_sum += rnd.Exp(1.);
       			L_sum += rnd.Lorentz(0., 1.);
       		}
       		write1 << (double) (U_sum / Ns[j]) << stampa[j];
       		write2 << E_sum / Ns[j] << stampa[j];
       		write3 << L_sum / Ns[j] << stampa[j];
     	}
	}	

	rnd.SaveSeed();
   	return 0;
}
