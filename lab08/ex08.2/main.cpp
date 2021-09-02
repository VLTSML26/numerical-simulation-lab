/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_grid		50
#define N_blocks	100

#define mu_min		0.5
#define mu_max		1.0
#define sigma_min	0.5
#define sigma_max	1.0

#include "../../Random/Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

double psi_T(double, double, double);
double V(double);
double D2_psi(double, double, double);

int main(int argc, char *argv[]) {

	// initialize random and output
	Random rnd;
	ofstream write;

	// define parameters, initialize x
	double mu, sigma;
	double x = 1.5;

	// delta has already been tuned: note that
	// you require to define Psi_T2 in order to tune it again.
	// This is different from lab08/ex08.1, since now Psi_T2 is
	// defined every time we iterate over the grid of (mu, sigma)
	double delta = 1.72;

	// runs the algorithm that finds the minimum of the energy
	if(string{argv[1]} == "Grid") {

		// less throws are required here
		const int N_throws = 1000;
		write.open("data/ex08.2_grid.out");
		mu = mu_min;
		sigma = sigma_min;
		double mu_step = (mu_max - mu_min) / N_grid;
		double sigma_step = (sigma_max - sigma_min) / N_grid;

		// 2-dimensional grid (mu, sigma)
		for(int i=0; i<N_grid; ++i) {
			mu = mu_min + i * mu_step;
			for(int j=0; j<N_grid; ++j) {
				sigma = sigma_min + j * sigma_step;

				// data blocking method for value of energy H
				for(int k=0; k<N_blocks; ++k) {

					//defines function to sample: |psi|^2
					auto psi_T2 = [mu, sigma](double x[]) {
						double psi = psi_T(*x, mu, sigma);
						return psi * psi;
					};
					// here I thought about tuning every time, in order to
					// maintain a perfect acceptance rate: it is too expensive
					double H = 0.;
					for(int l=0; l<N_throws; ++l) {
						rnd.Metropolis<double>(&x, 1, delta, psi_T2, "Uniform");
						H += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
					}
					H /= N_throws;
					write << H << "\t";
				}
				write << endl;
			}
		}
	}

	// runs the algorithm that, given the duple (mu, sigma) that
	// minimizes the energy, calculates it
	else if(string{argv[1]} == "Sample") {

		// a larger number of throws is needed here
		const int N_throws = 10000;
		write.open("data/ex08.2_measure.out");
		ofstream sample("data/ex08.2_sample.out");

		// reads mu and sigma from input. NOTE that I require to
		// run the code twice as I carried out the search for the 
		// minima of H with python in the notebook
		mu = atof(argv[2]);
		sigma = atof(argv[3]);

		// defines the function to be sampled
		auto psi_T2 = [mu, sigma](double x[]) {
			double psi = psi_T(*x, mu, sigma);
			return psi * psi;
		};

		// data blocking for H and also sampling of x from the distribution
		for(int i=0; i<N_blocks; ++i) {
			double H = 0.;
			for(int j=0; j<N_throws; ++j) {
				rnd.Metropolis<double>(&x, 1, delta, psi_T2, "Uniform");
				H += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
				sample << x << endl;
			}
			H /= N_throws;
			write << H << endl;
		}
	}

	// gives instruction if wrong input
	else {
		cerr << "Usage: " << argv[0] << " Grid >> if you want to find the minimum of H\n" <<
		"Otherwise\nUsage: " << argv[0] << " Sample <mu> <sigma> >> if you already know it\n";
		return -1;
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
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	phi_p *= phi_p;
	phi_m *= phi_m;
	double y = (phi_m - 1.) * exp(-0.5 * phi_m) + (phi_p - 1.) * exp(-0.5 * phi_p);
	return y / (sigma * sigma);
}
