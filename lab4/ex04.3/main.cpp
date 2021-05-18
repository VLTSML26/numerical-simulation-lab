/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "NVE.h"
#include "../../MolDyn/MolDyn.h"

using namespace std;

int main(int argc, char* argv[]) {

	int iprint = 100000;
	int imove = 1000;
	bool equilibrate, block;
	string old;
  	const char stampa[4] = {'\t', '\t', '\t', '\n'};
	
	// checking correct number of inputs from ./
  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}

	// choose start/equilibrate/measure
	if(string(argv[1]) == "start") {
		equilibrate = false;
		block = false;
	}
	else if(string(argv[1]) == "equilibrate") {
		equilibrate = true;
		block = false;
		old = "../input/old.0";
	}
	else if(string(argv[1]) == "measure") {
		equilibrate = false;
		block = true;
		old = "../input/old.0";
	}
	else {
		cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}
	
	// start simulation 
	NVE sim(equilibrate, old);
	ofstream write_blocks("block_measure.out", ios::app);
	
	// iprint = 10^3 nsteps = 10^3
  	for(int i=1; i <= sim.Get_blocks(); ++i) {
		double sum[sim.Get_props()] = {0.};
		for(int j=1; j <= sim.Get_throws(); ++j) {
			// cute counter
			if((i*sim.Get_throws() + j) % iprint == 0) 
				cout << (i*sim.Get_throws() + j)/iprint - 1 << "%\r";
			cout.flush();
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
				//sim.Measure();
			}
		}
		if(block) {
			for(int k=0; k<sim.Get_props(); ++k) {
				sum[k] /= sim.Get_throws();
				write_blocks << sum[k] << stampa[k];
			}
		}
  	}
  	
	write_blocks.close();
	cout << endl; 
	sim.ConfFinal("../input/config.final", "../input/old.final");

  	return 0;
}
