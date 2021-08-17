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

int main() {

	NVT sim;
	int iprint = 1000;
	char stampa[sim.Get_props()];
	for(int k=0; k<sim.Get_props() - 1; ++k) 
		stampa[k] = '\t';
	stampa[sim.Get_props()-1] = '\n';

	ofstream write("measure.out", ios::app);
	for(int i=1; i <= sim.Get_blocks(); ++i) {
		double sum[sim.Get_props()] = {0.};
		for(int j=1; j <= sim.Get_throws(); ++j) {
			// cute counter
			if((i*sim.Get_throws() + j) % iprint == 0) 
				cout << (i*sim.Get_throws() + j)/iprint - 1 << "%\r";
			cout.flush();
			sim.Move();
			sim.Measure();
			for(int k=0; k<sim.Get_props(); ++k)
				sum[k] += sim.walker(k);
		}
		sim.Measure(sum);
		for(int k=0; k<sim.Get_props(); ++k)
			write << sum[k] << stampa[k];
	}

	write.close();
	sim.ConfFinal();
	
	return 0;
}
