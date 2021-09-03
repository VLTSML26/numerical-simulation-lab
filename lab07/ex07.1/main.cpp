/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_equilibrate	500000
#define N_steps			10000
#define N_tune			1000
#define iprint			5000

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "NVT.h"

using namespace std;

int main() {

	// initialize simulation
	NVT sim;

	// tune metropolis
	cout << "==== Tuning Metropolis (acceptance rate should be near 0.6 due to PBC conditions)..." << endl;
	sim.Tune(N_tune);

	// equilibration
	cout << "==== Equilibration..." << endl;
	for(int i=0; i<N_equilibrate; ++i) {
		sim.Move();
	}

	// start measuring
	// NOTE that i sometimes use a different file
	// named ex07.1_big.out that has been produced
	// for N_steps = 500000
	cout << "==== Measure..." << endl;
	ofstream write("data/ex07.1.out", ios::app);
	for(int i=1; i<=N_steps; ++i) {
		sim.Move();
		sim.Measure();

		// write instantaneous values
		write << sim.walker(0) << "\t" << sim.walker(1) << endl;
	}

	cout << "Final Metropolis acceptance rate: " << (double) sim.Get_accepted() / N_steps / sim.Get_npart()<< endl;
	write.close();
	sim.ConfFinal();
	
	return 0;
}
