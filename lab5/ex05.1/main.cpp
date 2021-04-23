/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

template<class T>double distance(T*);
double psi_100(double[]);
double psi_210(double[]);

int main() {

	Random rnd;
	ofstream write_blocks, write_points;

  	const int N_throws = 100000;
  	const int N_blocks = 100;

	string sampling[2] = {"Uniform", "Gauss"};

	double x[2][3];
	for(int m=0; m<2; m++) {
		for(int n=0; n<3; n++) {
			x[m][n] = 1.5;
		}
	}
	double delta[2] = {1.2, 0.7};
	for(int m=0; m<2; m++)
		rnd.Tune(x[m], delta[m], psi_100, sampling[m]);

	write_blocks.open("ex05.1.out");
	write_points.open("ex05.1_points.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum[2] = {0.};
		for(int j=0; j<N_throws; j++) {
			for(int m=0; m<2; m++) {
				sum[m] += sqrt(distance<double>(x[m]));
				rnd.Metropolis(x[m], delta[m], psi_100, sampling[m]);
			//write_points << x[0] << "\t" << x[1] << "\t" << x[2] << endl; 
			}
		}
		for(int m=0; m<2; m++)
			sum[m] /= N_throws;
		write_blocks << sum[0] << "\t" << sum[1] << endl;
	}
	write_blocks.close();
	write_points.close();

	for(int m=0; m<2; m++) {
		for(int n=0; n<3; n++) {
			x[m][n] = 5.;
		}
	}
	delta[0] = 2.9;
	delta[1] = 1.8;
	for(int m=0; m<2; m++)
		rnd.Tune(x[m], delta[m], psi_210, sampling[m]);

	write_blocks.open("ex05.2.out");
	write_points.open("ex05.2_points.out");
   	for(int i=0; i<N_blocks; i++) {
		double sum[2] = {0.};
		for(int j=0; j<N_throws; j++) {
			for(int m=0; m<2; m++) {
				sum[m] += sqrt(distance<double>(x[m]));
				rnd.Metropolis(x[m], delta[m], psi_210, sampling[m]);
			//write_points << x[0] << "\t" << x[1] << "\t" << x[2] << endl; 
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

double psi_100(double x[]) {
	double r = sqrt(distance<double>(x));
	return exp(- 2. * r);
}

double psi_210(double x[]) {
	double r = sqrt(distance<double>(x));
	double psi = x[2] * exp(-r / 2.);
	return psi * psi;
}

template<class T> double distance(T* v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
