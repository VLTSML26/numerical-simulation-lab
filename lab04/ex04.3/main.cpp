/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define iprint	10000
#define imove	100

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "NVE.h"
#include "../../MolDyn/MolDyn.h"

using namespace std;

int main(int argc, char* argv[]) {

	// check correct number of inputs from terminal
  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}

	// prepare simulation for start/equilibrate/measure
	bool equilibrate, block;
	string old;
	string folder = "../input";
	if(string(argv[1]) == "start") {
		equilibrate = false;
		block = false;
	}
	else if(string(argv[1]) == "equilibrate") {
		equilibrate = true;
		block = false;
		old = folder + "/old.0";
	}
	else if(string(argv[1]) == "measure") {
		equilibrate = false;
		block = true;
		old = folder + "/old.0";
	}
	else {
		cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}
	
	// start simulation and initiaize output 
	NVE sim(folder, equilibrate, old);
	ofstream write_blocks("block_measure.out", ios::app);

	// trick in order to get the job done in one for cycle
  	char stampa[sim.Get_props()];
	for(int k=0; k<sim.Get_props() - 1; ++k) 
		stampa[k] = '\t';
	stampa[sim.Get_props()-1] = '\n';
	
  	for(int i=1; i <= sim.Get_blocks(); ++i) {
		double sum[sim.Get_props()] = {0.};
		for(int j=1; j <= sim.Get_throws(); ++j) {

			// cute counter (I leave it commented, it's not necessary)
			/*if((i*sim.Get_throws() + j) % iprint == 0) 
			 *	cout << (i*sim.Get_throws() + j)/iprint - 1 << "%\r";
			 *cout.flush();
			 */

			// during measure phase (block average)
			if(block) {
				sim.Move();
				sim.Measure();
				for(int k=0; k<sim.Get_props(); ++k)
					sum[k] += sim.walker(k);
			}

			// during start or equilibration phases one has to perform normal 
			// measures (not blk) in order check if the system has a good behavior
			// here commented in order to improve performance
			else if((i*sim.Get_throws() + j + 1) % imove == 0) {
				sim.Move();
//				sim.Measure();
			}
		}
		if(block) {

			// the first 4 elements of sum are the observables, need to divide by nthrows
			// no need to divide the remaining m_props - 4, used for g(r)
			for(int k=0; k<4; ++k)
				sum[k] /= (double) (sim.Get_throws() * sim.Get_npart());

			// print results
			for(int k=0; k<sim.Get_props(); ++k)
				write_blocks << sum[k] << stampa[k];
		}
  	}
  	
	write_blocks.close();
	cout << endl; 

	// print final configuration in order to restart
	sim.ConfFinal(folder + "/config.final", folder + "/old.final");

  	return 0;
}
