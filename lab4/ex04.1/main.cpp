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

using namespace std;

int main(int argc, char* argv[]) {

	int nstep = 10000;
  	const char stampa[4] = {'\t', '\t', '\t', '\n'};
	bool equilibrate;
	string old;

	// checking correct number of inputs from ./
  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
	  	return -1;
	}

	// choose start/equilibrate/measure
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
	
	// start simulation 
	NVE sim(equilibrate, old);
	ofstream write("measure.out", ios::app);

  	for(int i=1; i <= nstep; ++i) {
		sim.Move();
     	if(i%10 == 0) {
        	sim.Measure();
			for(int k=0; k<4; ++k)
				write << sim.walker(k) << stampa[k];
     	}
		// cute counter
		if(i%100 == 0) {
			cout << i/100 << "%\r";
			cout.flush();
		}
  	}
  	
	cout << endl;
	sim.ConfFinal("../input/config.final", "../input/old.final");
  	return 0;
}
