/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "NVT.h"

using namespace std;

int main(int argc, char *argv[]) {

	// initialize simulation
	string folder = "../input_NVT";
	NVT sim(folder);

	// read phase from input and initialize output
	string phase = argv[1];
	ofstream write("data/NVT_" + phase + ".out");

	// trick in order to produce only one output file
	char stampa[sim.Get_props()];
	for(int k=0; k<sim.Get_props() - 1; ++k) 
		stampa[k] = '\t';
	stampa[sim.Get_props()-1] = '\n';

	// data blocking
	for(int i=1; i<=sim.Get_blocks(); ++i) {
		double sum[sim.Get_props()] = {0.};
		for(int j=1; j<=sim.Get_throws(); ++j) {
			sim.Move();
			sim.Measure();

			// accumulates measures in sum
			for(int k=0; k<sim.Get_props(); ++k)
				sum[k] += sim.walker(k);
		}

		// tail corrections and rescaling bin occupations for g(r)
		sim.Measure(sum);

		// write results
		for(int k=0; k<sim.Get_props(); ++k)
			write << sum[k] << stampa[k];
	}

	write.close();
	sim.ConfFinal();
	
	return 0;
}
