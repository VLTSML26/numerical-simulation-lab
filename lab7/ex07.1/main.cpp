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

	int N_equilibrate = 500000;
	int N_steps = 10000;
	int N_tune = 1000;
	int iprint = (int) N_steps / 100;

	NVT sim;
	cout << "==== Tuning Metropolis (acceptance rate should be near 0.6 due to PBC conditions)..." << endl;
	sim.Tune(N_tune);
	cout << "==== Equilibration..." << endl;
	for(int i=0; i<N_equilibrate; ++i) {
		sim.Move();
///* cute counter
		if(i % (N_equilibrate / 100) == 0) 
			cout << i * 100 / N_equilibrate << "%\r";
		cout.flush();
//*/
	}

	cout << "==== Measure..." << endl;
	ofstream write("ex07.1.out", ios::app);
	for(int i=1; i <= N_steps; ++i) {
///* cute counter
		if(i % iprint == 0) 
			cout << i / iprint << "%\r";
		cout.flush();
//*/
		sim.Move();
		sim.Measure();
		write << sim.walker(0) << "\t" << sim.walker(1) << endl;
	}

	cout << "Final Metropolis acceptance rate: " << (double) sim.Get_accepted() / N_steps / sim.Get_npart()<< endl;
	write.close();
	sim.ConfFinal();
	
	return 0;
}
