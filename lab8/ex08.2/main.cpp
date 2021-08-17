/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "../../Random/Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

double psi_T(double, double, double);
double V(double);
double D2_psi(double, double, double);
/*double Dmu_psi(double, double, double);
double Dsigma_psi(double, double, double);
double Dmu_D2_psi(double, double, double);
double Dsigma_D2_psi(double, double, double);
double Dmu_H(double, double, double);
double Dsigma_H(double, double, double);
*/

int main() {

	Random rnd;
	ofstream write;

	const double mu_min = 1.;
	const double mu_max = 5.;
	const double sigma_min = 1.;
	const double sigma_max = 5.;
	const int N_grid = 300;
	const int N_steps = 10000;
  	const int N_throws = 100000;
  	const int N_blocks = 100;

	cout << "ciao" << endl;
	double delta = 3.525; // tuned with Random::Tune with values of mu_min and sigma_min
	double x = 1.5;

	double mu = mu_min;
	double sigma = sigma_min;
	auto psi_T2 = [mu, sigma](double x[]) {
		double psi = psi_T(x[0], mu, sigma);
		return psi * psi;
	};
	rnd.Tune(&x, 1, delta, psi_T2, "Uniform");
	cout << delta << endl;

	cout << "ciao" << endl;
	write.open("grid.out");
	for(int i=0; i<N_grid; ++i) {
		mu += (mu_max - mu_min) / N_grid;
		sigma = sigma_min;
		cout << i << endl;
		for(int j=0; j<N_grid; ++j) {
			sigma += (sigma_max - sigma_min) / N_grid;
//			rnd.Tune(&x, 1, delta, psi_T2, "Uniform");
			double H = 0.;
			for(int k=0; k<N_steps; ++k) {
				rnd.Metropolis(&x, 1, delta, psi_T2, "Uniform");
				H += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
			}
			H /= N_steps;
			write << mu << "\t" << sigma << "\t" << H << endl;
		}
	}
	write.close();

	/*
	// gradient descent
	int i = 0;
	double Dmu, Dsigma;
	do {
		Dmu = 0.;
	   	Dsigma = 0.;
		for(int j=0; j<N_steps; ++j) {
			rnd.Metropolis(&x, 1, delta, psi_T2, "Uniform");
			Dmu += Dmu_H(x, mu, sigma);
			Dsigma += Dsigma_H(x, mu, sigma);
		}
		Dmu /= N_steps;
		Dsigma /= N_steps;
		mu += Dmu * rate;
		sigma += Dsigma * rate;
		cout << i << "\t" << Dmu << "\t" << Dsigma << endl;
		i++;
		cout << (i<N_max) << endl;
		cout << (Dmu > 1e-3 or Dsigma > 1e-3) << endl;
	} while(i<N_max and (Dmu > 0.001 or Dsigma > 0.001));
*/
	/*
	write_blocks.open("ex08.2.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum = 0.;
		for(int j=0; j<N_throws; j++) {
			sum += (-0.5 * D2_psi(x, mu, sigma)) / psi_T(x, mu, sigma) + V(x);
			rnd.Metropolis(&x, 1, delta, psi_T2, "Uniform");
		}
		sum /= N_throws;
		write_blocks << sum << endl;
	}
	write_blocks.close();
*/
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
/*
double Dmu_psi(double x, double mu, double sigma) {
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	double psi_p = exp(-0.5 * phi_p * phi_p);
	double psi_m = exp(-0.5 * phi_m * phi_m);
	return (phi_m*psi_m - phi_p*psi_p) / sigma;
}

double Dsigma_psi(double x, double mu, double sigma) {
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	double psi_p = exp(-0.5 * phi_p * phi_p);
	double psi_m = exp(-0.5 * phi_m * phi_m);
	return (phi_m*phi_m*psi_m + phi_p*phi_p*psi_p) / sigma;
}

double Dmu_D2_psi(double x, double mu, double sigma) {
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	double psi_p = exp(-0.5 * phi_p * phi_p);
	double psi_m = exp(-0.5 * phi_m * phi_m);
	return ((-1. - phi_m + phi_m*phi_m)*psi_m + (1. + phi_p - phi_p*phi_p)*psi_p) / sigma;
}

double Dsigma_D2_psi(double x, double mu, double sigma) {
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	double psi_p = exp(-0.5 * phi_p * phi_p);
	double psi_m = exp(-0.5 * phi_m * phi_m);
	return ((1. - phi_m + phi_m*phi_m)*psi_m*phi_m + (1. - phi_p + phi_p*phi_p)*psi_p*phi_p) / sigma;
}

double Dmu_H(double x, double mu, double sigma) {
	return Dmu_D2_psi(x, mu, sigma) / psi_T(x, mu, sigma) - (D2_psi(x, mu, sigma)*Dmu_psi(x, mu, sigma)) / (psi_T(x, mu, sigma)*psi_T(x, mu, sigma));
}

double Dsigma_H(double x, double mu, double sigma) {
	return Dsigma_D2_psi(x, mu, sigma) / psi_T(x, mu, sigma) - (D2_psi(x, mu, sigma)*Dsigma_psi(x, mu, sigma)) / (psi_T(x, mu, sigma)*psi_T(x, mu, sigma));
}*/
