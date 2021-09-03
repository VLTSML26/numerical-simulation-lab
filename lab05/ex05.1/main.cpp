/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_throws	100000
#define N_blocks	100

#include "../../Random/Random.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

template<class T>double distance(T*);
double psi_100(double[]);
double psi_210(double[]);

int main(int argc, char* argv[]) {

	// initialize random and outputs
	Random rnd;
	ofstream write_blocks, write_points;

	string sampling[2] = {"Uniform", "Gauss"};

	// start simulation for psi_100
	cout << "====== PSI_100 ======" << endl;

	// read distance from origin (and if it is
	// considered near of far) from terminal
	double r = atof(argv[1]);
	string type = argv[3];

	// randomly place 2 points (x,y,z) in space at fixed r
	double x[2][3];
	for(int m=0; m<2; m++)
		rnd.Spherical3D(r, x[m]);

	// tune metropolis in order to have acceptance rate 0.5
	// (it should be already tuned)
	double delta[2] = {1.22, 0.76};
	for(int m=0; m<2; m++) {
		rnd.Tune<double>(x[m], 3, delta[m], psi_100, sampling[m]);
		cout << "Tuned " << sampling[m] << " with value delta = " << delta[m] << endl;
	}

	// data blocking
	write_blocks.open("data/" + type + "_psi100.out");
	write_points.open("data/" + type + "_psi100_points.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum[2] = {0.};
		for(int j=0; j<N_throws; j++) {
			for(int m=0; m<2; m++) {
				sum[m] += sqrt(distance<double>(x[m]));
				rnd.Metropolis<double>(x[m], 3, delta[m], psi_100, sampling[m]);

				// prints the points in order to visualize the orbital
				if(j % 1000 == 0)
					write_points << x[m][0] << "\t" << x[m][1] << "\t" << x[m][2] << endl; 
			}
		}
		for(int m=0; m<2; m++)
			sum[m] /= N_throws;
		write_blocks << sum[0] << "\t" << sum[1] << endl;
	}
	write_blocks.close();
	write_points.close();

	// start simulation for psi_210
	cout << "====== PSI_210 ======" << endl;

	// read distance from origin from terminal
	r = atof(argv[2]);

	for(int m=0; m<2; m++)
		rnd.Spherical3D(r, x[m]);

	// tune metropolis in order to have acceptance rate 0.5
	// (it should be already tuned)
	delta[0] = 2.97;
	delta[1] = 1.87;
	for(int m=0; m<2; m++) {
		rnd.Tune<double>(x[m], 3, delta[m], psi_210, sampling[m]);
		cout << "Tuned " << sampling[m] << " with value delta = " << delta[m] << endl;
	}

	// data blocking
	write_blocks.open("data/" + type + "_psi210.out");
	write_points.open("data/" + type + "_psi210_points.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum[2] = {0.};
		for(int j=0; j<N_throws; j++) {
			for(int m=0; m<2; m++) {
				sum[m] += sqrt(distance<double>(x[m]));
				rnd.Metropolis<double>(x[m], 3, delta[m], psi_210, sampling[m]);

				// prints the points in order to visualize the orbital
				if(j % 100 == 0)
					write_points << x[m][0] << "\t" << x[m][1] << "\t" << x[m][2] << endl; 
			}
		}
		for(int m=0; m<2; m++)
			sum[m] /= N_throws;
		write_blocks << sum[0] << "\t" << sum[1] << endl;
	}
	write_blocks.close();
	write_points.close();
	rnd.SaveSeed();
   	return 0;
}

/*
 * This is the GS wavefunction psi_100
 */ 
double psi_100(double x[]) {
	double r = sqrt(distance<double>(x));
	return exp(- 2. * r);
}

/*
 * This is the excited wavefunction psi_210
 */
double psi_210(double x[]) {
	double r = sqrt(distance<double>(x));
	double psi = x[2] * exp(-r / 2.);
	return psi * psi;
}

/*
 * This function returns the L2 distance from the origin
 */
template<class T> double distance(T* v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
