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
#include <string>
#include <cmath>

using namespace std;

double psi_T(double, double, double);
double V(double);
double D2_psi(double, double, double);

int main() {

	// initialize random and output
	Random rnd;
	ofstream write_blocks("data/ex08.1.out");

	// initialize delta (already tuned), x and parameters
	double delta = 2.365;
	double x = 1.5;
	const double mu = 0.5, sigma = 1.;

	// defines the function to sample: |psi|^2
	auto psi_T2 = [mu, sigma](double x[]) {
		double psi = psi_T(x[0], mu, sigma);
		return psi * psi;
	};

	// this should take only few steps since delta is already tuned
	rnd.Tune<double>(&x, 1, delta, psi_T2, "Uniform");

	// data blocking for the energy
   	for(int i=0; i<N_blocks; i++) {
		double sum = 0.;
		for(int j=0; j<N_throws; j++) {
			// energy
			sum += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
			// samples new x
			rnd.Metropolis<double>(&x, 1, delta, psi_T2, "Uniform");
		}
		sum /= N_throws;
		write_blocks << sum << endl;
	}

	rnd.SaveSeed();
   	return 0;
}

/*
 * This is the trial wavefunction psi_T(x)
 */
double psi_T(double x, double mu, double sigma) {
	double phi_p = x + mu;
	double phi_m = x - mu;
	double y = 0.5 / (sigma * sigma);
	phi_p *= phi_p;
	phi_m *= phi_m;
	return exp(-y * phi_p) + exp(-y * phi_m);
}

/*
 * This is the potential V(x)
 */
double V(double x) {
	return x*x*x*x - 2.5 * x*x;
}

/*
 * This is the second derivative of the trial wavefunction
 */
double D2_psi(double x, double mu, double sigma) {
	double phi_p = x + mu;
	double phi_m = x - mu;
	double y = 0.5 / (sigma * sigma);
	phi_p *= phi_p;
	phi_m *= phi_m;
	return (phi_m - 1.) * exp(-y * phi_m) + (phi_p - 1.) * exp(-y * phi_p);
}
