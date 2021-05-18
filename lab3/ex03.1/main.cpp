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
  	const int N_path = 100; // per dividere l'intervallo [0, t1] in N_path bin di egual lunghezza (ex03.2)
	const double s0 = 100.; // asset price at t=0
	const double k = 100.; // strike price
	const double sigma = 0.25; // volatility
	const double r = 0.1; // risk-free interest rate
	const double t1 = 1.; // delivery time

	// ==== EX03.1: DIRECT SAMPLE ====
   	
	ofstream write1("ex03.1.1.out");
	for(int i=0; i<N_blocks; i++) {
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; j++) {
			double z = rnd.Gauss(0., 1.);
			double s = s0 * exp((r - sigma*sigma/2.)*t1 + sigma*z*sqrt(t1));
			s -= k;
			if(s > 0) 
				C_sum += s*exp(- r*t1); // il profitto per una Call è max(0, s - k) ovvero max(0, s)
			else 
				P_sum += - s*exp(- r*t1); // il profitto per una Put è max(0, k - s) ovvero max(0, -s)
		}
		C_sum /= N_throws;
		P_sum /= N_throws;
     	write1 << C_sum << "\t" << P_sum << endl;
	}
	write1.close();

	// ==== EX03.2: DISCRETIZED SAMPLE ====

	ofstream write2("ex03.1.2.out");
	for(int i=0; i<N_blocks; i++) {
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; j++) {
			double s = s0;
			for(int k=0; k<N_path; k++) {
				double z = rnd.Gauss(0., 1.);
				s *= exp((r - sigma*sigma/2.)*t1/N_path + sigma*z*sqrt(t1/N_path));
			}
			s -= k;
			if(s > 0) 
				C_sum += s*exp(- r*t1);
			else 
				P_sum += - s*exp(- r*t1);
		}
		C_sum /= N_throws;
		P_sum /= N_throws;
     	write2 << C_sum << "\t" << P_sum << endl;
	}
	write2.close();
   	
	rnd.SaveSeed();
   	return 0;
}
