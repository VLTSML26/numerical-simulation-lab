/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_RWs		1000
#define N_steps		100
#define N_blocks	100

#include "../../Random/Random.h"
#include <fstream>
#include <cmath>

using namespace std;
 
template<class T>double distance(T*);
void lattice_RW(int*, double*, Random&);
void continuum_RW(double*, double*, Random&);

int main() {

	// initialize random and outputs
	Random rnd;
  	ofstream lattice ("data/ex02.2_lattice.out");
  	ofstream continuum ("data/ex02.2_continuum.out");

	for(int i=0; i<N_blocks; ++i) {

		// stores the distance at each step
		double L_r2[N_steps] = {0.};
		double C_r2[N_steps] = {0.};

     	for(int j=0; j<N_RWs; ++j) {
     		int L_r[3] = {0};
     		double C_r[3] = {0.};

			// calls algorithms for discrete and continuous RWs
     		lattice_RW(L_r, L_r2, rnd);
			continuum_RW(C_r, C_r2, rnd);
		}

		// writes data
		for(int k=0; k<N_steps; ++k) {
			L_r2[k] /= N_RWs;
			C_r2[k] /= N_RWs;

			// each line contains the distance at each step
			// the final value is what you need in order to
			// use data blocking and calculate <R>
			lattice << sqrt(L_r2[k]) << "\t";
			continuum << sqrt(C_r2[k]) << "\t";
		}
		lattice << endl;
		continuum << endl;
	}

	rnd.SaveSeed();
   	return 0;
}

/*
 * This function performs RW on a lattice
 */
void lattice_RW(int* r, double* r2, Random &rnd) {
   	for(int k=0; k<N_steps; ++k) {

		// randomly chooses a direction on the lattice
		// 0 -> x
		// 1 -> y
		// 2 -> z
		int ran = (int) rnd.Rannyu(0, 3);

		// randomly chooses to take a step back of forward
		// true -> forward
		// false -> back
   		bool ran2 = rnd.Bool();

		// moves on the lattice
		if(ran == 0) {
       		if(ran2) 
				r[0]++;
       		else 
				r[0]--;
       	}
       	if(ran == 1) {
       		if(ran2) 
				r[1]++;
       		else 
				r[1]--;
       	}
    	if(ran == 2) {
       		if(ran2) 
				r[2]++;
       		else 
				r[2]--;
       	}

		// stores the distance each step
		r2[k] += distance<int>(r);
    }
}

/*
 * This function performs continuous RW in 3D space
 */
void continuum_RW(double* r, double *r2, Random &rnd) {
	for(int k=0; k<N_steps; ++k) {

		// randomly chooses a direction sampling theta and phi
		double phi = rnd.Rannyu(0, 2*M_PI);
		double theta = rnd.Rannyu(0, M_PI);

		// takes a step in spherical coordinates (r=1)
		r[0] += sin(theta) * cos(phi);
		r[1] += sin(theta) * sin(phi);
		r[2] += cos(theta);

		// stores the distance each step
		r2[k] += distance<double>(r);
	}
}

/*
 * This function returns the distance of a point from the origin
 */
template<class T> double distance(T* v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
