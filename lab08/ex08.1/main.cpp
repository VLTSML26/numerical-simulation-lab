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

	Random rnd;
	ofstream write_blocks;

	double mu = 0.5, sigma = 1.;
	auto psi_T2 = [mu, sigma](double x[]) {
		double psi = psi_T(x[0], mu, sigma);
		return psi * psi;
	};

	double delta = 2.365; // tuned with Random::Tune
	double x = 1.5;
	rnd.Tune(&x, 1, delta, psi_T2, "Uniform");

	write_blocks.open("ex08.1bis.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum = 0.;
		for(int j=0; j<N_throws; j++) {
//			sum += (-0.5 * D2_psi(x) + V(x) * psi_T(x)) / psi_T(x);
			sum += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
			rnd.Metropolis(&x, 1, delta, psi_T2, "Uniform");
		}
		sum /= N_throws;
		write_blocks << sum << endl;
	}
	write_blocks.close();

	rnd.SaveSeed();
   	return 0;
}

double psi_T(double x, double mu, double sigma) {
	double phi_p = x + mu;
	double phi_m = x - mu;
	double y = 0.5 / (sigma * sigma);
	phi_p *= phi_p;
	phi_m *= phi_m;
	return exp(-y * phi_p) + exp(-y * phi_m);
}

double V(double x) {
	return x*x*x*x - 2.5 * x*x;
}

double D2_psi(double x, double mu, double sigma) {
	double phi_p = x + mu;
	double phi_m = x - mu;
	double y = 0.5 / (sigma * sigma);
	phi_p *= phi_p;
	phi_m *= phi_m;
	return (phi_m - 1.) * exp(-y * phi_m) + (phi_p - 1.) * exp(-y * phi_p);
}
