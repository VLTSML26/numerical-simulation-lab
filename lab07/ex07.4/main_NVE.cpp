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

	// checks correct number of inputs from terminal and sets phase
  	if(argc<3) {
	  	cerr << "Usage: " << argv[0] << " <start/equilibrate/measure> <solid/liquid/gas>\n";
		return 1;
	}
	string phase = argv[2];

	// prepare simulation for start/equilibrate/measure
	bool equilibrate, block;
	string old;
	string folder = "../input_NVE";
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
	
	// start simulation 
	NVE sim(folder, equilibrate, old);
	ofstream write_blocks("data/NVE_" + phase + ".out", ios::app);
  	char stampa[sim.Get_props()];
	for(int k=0; k<sim.Get_props() - 1; ++k) 
		stampa[k] = '\t';
	stampa[sim.Get_props()-1] = '\n';
	
  	for(int i=1; i <= sim.Get_blocks(); ++i) {
		double sum[sim.Get_props()] = {0.};
		for(int j=1; j <= sim.Get_throws(); ++j) {

			// cute counter
			/*
			if((i*sim.Get_throws() + j) % iprint == 0) 
				cout << (i*sim.Get_throws() + j)/iprint - 1 << "%\r";
			cout.flush();
			*/

			// during measure phase (block average)
			if(block) {
				sim.Move();
				sim.Measure();
				for(int k=0; k<sim.Get_props(); ++k)
					sum[k] += sim.walker(k);
			}

			// during start or equilibration phases one has to perform normal measures (not blk)
			// in order check if the system has a good behavior
			// since it is already done in the previous exercises, I leave it commented
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

			// updates the g(r) estimation before writing output files (more in NVE.cpp)
			sim.Gofr(sum);

			// print results
			for(int k=0; k<sim.Get_props(); ++k)
				write_blocks << sum[k] << stampa[k];
		}
  	}
  	
	write_blocks.close();
	cout << endl; 

	// write final configuration
	sim.ConfFinal(folder + "/config.final", folder + "/old.final");

  	return 0;
}
