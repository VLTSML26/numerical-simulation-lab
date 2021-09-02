/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_throws	 10000
#define N_blocks	 100
#define N_path		 100	// in order to divide [0, t1] in N_path bins of equal lenght
#define s0			 100.	// asset price at t=0
#define k			 100.	// strike price
#define sigma		 0.25	// volatility
#define r			 0.1	// risk-free interest rate
#define t1			 1.		// delivery time

#include "../../Random/Random.h"
#include <fstream>
#include <cmath>

using namespace std;

int main() {

	// initialize random and outputs
	Random rnd;
	ofstream direct("data/ex03.1_direct.out");
	ofstream discretized("data/ex03.1_discretized.out");

	// data blocking method using direct sample 
	for(int i=0; i<N_blocks; ++i) {
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; ++j) {
			double z = rnd.Gauss(0., 1.);
			double s = s0 * exp((r - sigma*sigma/2.)*t1 + sigma*z*sqrt(t1));
			s -= k;
			if(s > 0) 

				// profit for a call is max(0, s-k), hence max(0, s)
				C_sum += s*exp(- r*t1);
			else 

				// profit for a put is max(0, k-s), hence max(0, -s)
				P_sum += - s*exp(- r*t1);
		}
		C_sum /= N_throws;
		P_sum /= N_throws;
     	direct << C_sum << "\t" << P_sum << endl;
	}

	// data blocking method using discretized sample 
	for(int i=0; i<N_blocks; ++i) {
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; ++j) {
			double s = s0;
			for(int n=0; n<N_path; ++n) {
				double z = rnd.Gauss(0., 1.);
				s *= exp((r - sigma*sigma/2.)*t1/N_path + sigma*z*sqrt(t1/N_path));
			}
			s -= k;
			if(s > 0) 

				// profit for a call is max(0, s-k), hence max(0, s)
				C_sum += s*exp(- r*t1);
			else 

				// profit for a put is max(0, k-s), hence max(0, -s)
				P_sum += - s*exp(- r*t1);
		}
		C_sum /= N_throws;
		P_sum /= N_throws;
     	discretized << C_sum << "\t" << P_sum << endl;
	}
   	
	rnd.SaveSeed();
   	return 0;
}
