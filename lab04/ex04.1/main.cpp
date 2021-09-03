/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_steps	10000

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "NVE.h"

using namespace std;

int main(int argc, char* argv[]) {

	// check correct number of inputs from terminal
  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
	  	return -1;
	}

	// prepare simulation for start/equilibrate/measure
	bool equilibrate;
	string old;
	if(string(argv[1]) == "start") 
		equilibrate = false;
	else if(string(argv[1]) == "equilibrate") {
		equilibrate = true;
		old = "old.0";
	}
	else if(string(argv[1]) == "measure") {
		equilibrate = false;
		old = "old.0";
	}
	else {
		cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}
	
	// start simulation and initialize output 
	NVE sim(equilibrate, old);
	ofstream write("data/measure.out", ios::app);

	// trick in order to get the job done in one for cycle
  	const char stampa[4] = {'\t', '\t', '\t', '\n'};

  	for(int i=1; i<=N_steps; ++i) {
		sim.Move();

		// it is not necessary to measure every time I move
		// the system, since N_steps is very big
     	if(i%10 == 0) {
        	sim.Measure();
			for(int k=0; k<4; ++k)
				write << sim.walker(k) << stampa[k];
     	}
		
		// cute counter (I leave it commented, it's not necessary)
		/*if(i%100 == 0) {
		 *	cout << i/100 << "%\r";
		 *	cout.flush();
		}*/
  	}
  	
	cout << endl;

	// print final configuration in order to restart
	sim.ConfFinal("../input/config.final", "../input/old.final");
  	return 0;
}
