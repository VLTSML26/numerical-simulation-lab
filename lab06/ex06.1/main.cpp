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
#include "Ising.h"

using namespace std;

int main(int argc, char* argv[]) {

	// checks correct number of inputs from terminal
  	if(argc < 4) {
	  	cerr << "Usage: " << argv[0] << " metropolis/gibbs <N_measures> <external field>\n";
	  	return -1;
	}

	// reads the sampling method to be used from terminal
	bool metro;
	if(string(argv[1]) == "metropolis") {
		metro = true;
	}
	else if(string(argv[1]) == "gibbs") {
		metro = false;
	}
	else {
	  	cerr << "Usage: " << argv[0] << " metropolis/gibbs <N_measure> <exernal field>\n";
	  	return 1;
	}

	// reads N_measure and magnetization h from terminal
	int N_measure = atoi(argv[2]);
	double h = atof(argv[3]);


	for(int i=0; i<N_measure; ++i) {

		// the following lines (commented) are just counters
		/*
		int percent = 100 * i / N_measure;
		cout << "=== " << percent << "%\r";
		cout.flush();
		*/

		// increment temperature and start simulation at fixed T
		double temp = 0.5 + 1.5 * (double) i / (N_measure - 1);
		Ising sim(metro, temp, h);

		// first step: need to equilibrate
		if(i == 0) 
			sim.Equilibrate();

		// data blocking
		for(int iblk=1; iblk<=sim.Get_blk(); ++iblk) {
			sim.Reset(iblk);
			for(int istep=1; istep <= sim.Get_step(); ++istep) {
				sim.Move();
				sim.Measure();
				sim.Accumulate();
			}
			sim.Averages(iblk);
		}

		// print final configuration and printf info
		sim.ConfFinal();
		if(i == N_measure - 1) {
			if(h == 0.0)
				cout << "Final results in *.out: 1-U, 2-C, 3-X, 4-M\n"; 
			else 
				cout << "Final results in *_extfield.out: 1-U, 2-C, 3-X, 4-M\n";
		}
	}

	return 0;
}
