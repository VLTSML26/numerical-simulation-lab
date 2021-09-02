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

	Random rnd;
	ofstream write;
	double mu, sigma;
	double x = 1.5;
	double delta = 1.72;	// tuned with Random::Tune<double>, now commented, with values of mu_min and sigma_min
							// note that you require to define Psi_T2 in order to tune again
	if(string{argv[1]} == "Grid") {
		const int N_throws = 1000;
		write.open("data/ex08.2_grid.out");
		mu = mu_min;
		sigma = sigma_min;
		double mu_step = (mu_max - mu_min) / N_grid;
		double sigma_step = (sigma_max - sigma_min) / N_grid;
		for(int i=0; i<N_grid; ++i) {
			mu = mu_min + i * mu_step;
			for(int j=0; j<N_grid; ++j) {
				sigma = sigma_min + j * sigma_step;
				for(int k=0; k<N_blocks; ++k) {
					auto psi_T2 = [mu, sigma](double x[]) {
						double psi = psi_T(*x, mu, sigma);
						return psi * psi;
					};
					// I thought about tuning every time, in order to maintain acceptance rate. It is simply too expensive
					//rnd.Tune<double>(&x, 1, delta, psi_T2, "Uniform");
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
	else if(string{argv[1]} == "Sample") {
		const int N_throws = 10000;
		write.open("data/ex08.2_measure.out");
		ofstream sample("data/ex08.2_sample.out");
		mu = atof(argv[2]);
		sigma = atof(argv[3]);
		cout << mu << "\t" << sigma << endl;
		auto psi_T2 = [mu, sigma](double x[]) {
			double psi = psi_T(*x, mu, sigma);
			return psi * psi;
		};
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
	else {
		cerr << "Usage: " << argv[0] << "<Grid/Sample>" << endl;
		return -1;
	}

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
	double phi_p = (x + mu) / sigma;
	double phi_m = (x - mu) / sigma;
	phi_p *= phi_p;
	phi_m *= phi_m;
	double y = (phi_m - 1.) * exp(-0.5 * phi_m) + (phi_p - 1.) * exp(-0.5 * phi_p);
	return y / (sigma * sigma);
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
