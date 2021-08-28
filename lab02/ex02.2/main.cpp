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
#define N_steps		10000
#define N_blocks	100

#include "../../Random/Random.h"
#include <fstream>
#include <cmath>

using namespace std;
 
template<class T>double distance(T*);
void lattice_RW(int*, Random&);
void continuum_RW(double*, Random&);

int main() {

	Random rnd;
  	ofstream write ("data/ex02.2.out");

	for(int i=0; i<N_blocks; ++i) {
     	double L_d = 0.;
     	double C_d = 0.;
     	for(int j=0; j<N_RWs; ++j) {
     		int L_r[3] = {0};
     		double C_r[3] = {0.};
     		lattice_RW(L_r, rnd); // compio il RW sul lattice: vedi sotto
			continuum_RW(C_r, rnd); // compio il RW nel continuo: vedi sotto
     		L_d += distance<int>(L_r);
     		C_d += distance<double>(C_r);
		}
		L_d /= N_RWs;
		C_d /= N_RWs;
     	write << L_d << "\t" << C_d << endl;
	}
	write.close();

	rnd.SaveSeed();
   	return 0;
}

void lattice_RW(int* r, Random &rnd) {
    	for(int k=0; k<N_steps; ++k) {
		int ran = (int) rnd.Rannyu(0, 3); // genero un numero tra {1,2,3} corrispondente a una direzione sul lattice
   		bool ran2 = rnd.Bool(); // true -> passo avanti, false -> passo indietro
		if(ran == 0) {
       			if(ran2) r[0]++;
       			else r[0]--;
       		}
       		if(ran == 1) {
       			if(ran2) r[1]++;
       			else r[1]--;
       		}
       		if(ran == 2) {
       			if(ran2) r[2]++;
       			else r[2]--;
       		}
     	}
}

void continuum_RW(double* r, Random &rnd) {
	for(int k=0; k<N_steps; ++k) {
		double phi = rnd.Rannyu(0, 2*M_PI);
		double theta = rnd.Rannyu(0, M_PI);
		r[0] += sin(theta) * cos(phi);
		r[1] += sin(theta) * sin(phi);
		r[2] += cos(theta);
	}
}

template<class T> double distance(T* v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
