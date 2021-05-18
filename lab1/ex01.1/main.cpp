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

  	ofstream write ("ex01.1.out");
	for(int i=0; i<N_blocks; i++) {
	    double sum[2] = {0.};
		int chi[N_blocks] = {0};
		for(int j=0; j<N_throws; j++) {
			double ran = rnd.Rannyu();
			sum[0] += ran;
			sum[1] += (ran - 0.5) * (ran - 0.5);
			chi[(int) (ran*N_blocks)] ++;
		}
		for(int k=0; k<2; k++)
			sum[k] /= N_throws;
		double chi2 = 0.;
		for(int k=0; k<N_blocks; k++) {
			chi2 += (chi[k] - N_throws/N_blocks) * (chi[k] - N_throws/N_blocks);
		}
		chi2 *= ((double) N_blocks / N_throws);
		write << sum[0] << "\t" << sum[1] << "\t" << chi2 << endl;
	}	

	rnd.SaveSeed();
   	return 0;
}
