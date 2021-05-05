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

  	if(argc < 3) {
	  	cerr << "Usage: " << argv[0] << " metropolis/gibbs <N_measures>\n";
	  	return -1;
	}

	bool metro;

	if(string(argv[1]) == "metropolis") {
		metro = true;
	}
	else if(string(argv[1]) == "gibbs") {
		metro = false;
	}
	else {
	  	cerr << "Usage: " << argv[0] << " metropolis/gibbs <N_measure>\n";
	  	return 1;
	}

	int N_measure = atoi(argv[2]);

	for(int i=0; i<N_measure; ++i) {
		int percent = 100 * i / N_measure;
		cout << "=== " << percent << "%\r";
		cout.flush();
		double temp = 0.5 + 1.5 * (double) i / (N_measure - 1);
		Ising sim(metro, temp);
		sim.Equilibrate();
		for(int iblk=1; iblk <= sim.Get_blk(); ++iblk) {
			sim.Reset(iblk);
			for(int istep=1; istep <= sim.Get_step(); ++istep) {
				sim.Move();
				sim.Measure();
				sim.Accumulate();
			}
			sim.Averages(iblk);
		}
		sim.ConfFinal();
		if(i == N_measure - 1) {
			if(sim.Get_h() == 0.0)
				cout << "Final results in *.out: 1-U, 2-C, 3-M, 4-X\n"; 
			else 
				cout << "Final results in *_extfield.out: 1-U, 2-C, 3-M, 4-X\n";
		}
	}


	return 0;
}
