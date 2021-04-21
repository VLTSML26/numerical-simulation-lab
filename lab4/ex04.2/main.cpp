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

	bool equilibrate, block;
	string old;

	if(string(argv[1]) == "start") {
		equilibrate = false;
		block = false;
	}
	else if(string(argv[1]) == "equilibrate") {
		equilibrate = true;
		block = false;
		old = "old.0";
	}
	else if(string(argv[1]) == "measure") {
		equilibrate = false;
		block = true;
		old = "old.0";
	}
	else {
		cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}
	
	MolDyn sim(equilibrate, old);
	ofstream write_blocks("block_measure.out", ios::app);
	
  	for(int i=1; i <= sim.Get_blocks(); ++i) {
		double sum[4] = {0.};
		for(int j=1; j <= sim.Get_throws(); ++j) {
			sim.Move();
			if((i*sim.Get_throws() + j) % sim.Get_iprint() == 0) 
				cout << (i*sim.Get_throws() + j)/sim.Get_iprint() * 10 << "%\r";
			cout.flush();
			if(block) {
				double appo[4] = {0.};
				sim.Measure(appo);
				for(int k=0; k<4; ++k)
					sum[k] += appo[k];
			}
			else if((i*sim.Get_throws() + j + 1) % 10 == 0)
				sim.Measure();
		}
		if(block) {
			for(int k=0; k<4; ++k)
				sum[k] /= sim.Get_throws();
			write_blocks << sum[0] << "\t" << sum[1] << "\t"
						 << sum[2] << "\t" << sum[3] << "\n";
		}
  	}
  	
	write_blocks.close();
	cout << endl;
	sim.ConfFinal("config.final", "old.final");
  	return 0;
}
