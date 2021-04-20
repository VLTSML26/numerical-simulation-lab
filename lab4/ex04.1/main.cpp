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
#include "MolDyn.h"

using namespace std;

int main(int argc, char* argv[]) {

  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
	  	return -1;
	}

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
	
	MolDyn sim(equilibrate, old);
//	int nconf = 1;
  	for(int i=1; i <= sim.Get_steps(); ++i) {
		sim.Move();
     	if(i % sim.Get_iprint() == 0) cout << i/sim.Get_iprint() * 10 << "%\r";
		cout.flush();
     	if(i%10 == 0) {
        	sim.Measure();
//	        ConfXYZ(nconf); 
//     		nconf ++;
     	}
  	}
  	
	cout << endl;
	sim.ConfFinal("config.final", "old.final");
  	return 0;
}
