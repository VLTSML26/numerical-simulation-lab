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

	// initialize random and output
	Random rnd;
	ofstream write("data/ex01.1.out");

	// data blocking method
	for(int i=0; i<N_blocks; ++i) {

		// prepare 2D array
		// 0 - integral <r>
		// 1 - sigma^2
	    double sum[2] = {0.};

		// N_blocks-dimensional array for chi^2
		int chi[N_blocks] = {0};

		for(int j=0; j<N_throws; ++j) {
			double ran = rnd.Rannyu();
			sum[0] += ran;
			sum[1] += (ran - 0.5) * (ran - 0.5);
			chi[(int) (ran*N_blocks)] ++;
		}
	
		// prepare data for output: here i only divide for
		// N_throws, while the data analysis is carried
		// out in python
		for(int k=0; k<2; ++k)
			sum[k] /= N_throws;
		double chi2 = 0.;
		for(int k=0; k<N_blocks; ++k) {
			chi2 += (chi[k] - N_throws/N_blocks) * (chi[k] - N_throws/N_blocks);
		}
		chi2 *= ((double) N_blocks / N_throws);
		write << sum[0] << "\t" << sum[1] << "\t" << chi2 << endl;
	}	

	rnd.SaveSeed();
   	return 0;
}
