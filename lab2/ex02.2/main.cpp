/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
 
double error(double, double, int);
template<class T>double distance(T*);
void lattice_RW(int*, int, Random&);
void continuum_RW(double*, int, Random&);
Random init();

int main (int argc, char *argv[]){

	Random rnd;
	rnd = init();

  	const int N_RWs = 1000;
  	const int N_steps = 10000;
  	const int N_blocks = 100;

  	ofstream data1 ("ex02.2.1.out"); // lattice RW
	ofstream data2 ("ex02.2.2.out"); // continuum RW

	// ==== NOTA BENE ====
	// L sta per Lattice RW
	// C sta per Continuum RW

	double L_sum_prog = 0., L_sum2_prog = 0.;
	double C_sum_prog = 0., C_sum2_prog = 0.;	
   	
	for(int i=0; i<N_blocks; i++){
     		double L_d = 0.;
     		double C_d = 0.;
     		for(int j=0; j<N_RWs; j++){
     			int L_r[3] = {0};
     			double C_r[3] = {0.};
     			lattice_RW(L_r, N_steps, rnd); // compio il RW sul lattice: vedi sotto
			continuum_RW(C_r, N_steps, rnd); // compio il RW nel continuo: vedi sotto
     			L_d += distance<int>(L_r);
     			C_d += distance<double>(C_r);
		}
		// ==== ANALISI DATI (metodo a blocchi) ====
		// L
		L_sum_prog += L_d / N_RWs;
     		L_sum2_prog += L_d / N_RWs * L_d / N_RWs;
     		double av = sqrt(L_sum_prog / (i + 1)); // N.B. che serve la radq stavolta
     		double av2 = sqrt(L_sum2_prog / (i + 1));
     		data1 << av << "\t" << error(av, av2, i) << endl;
     		// C
     		C_sum_prog += C_d / N_RWs;
     		C_sum2_prog += C_d / N_RWs * C_d / N_RWs;
     		av = sqrt(C_sum_prog / (i + 1));
     		av2 = sqrt(C_sum2_prog / (i + 1)); 
     		data2 << av << "\t" << error(av, av2, i) << endl;
	}

	rnd.SaveSeed();
   	return 0;
}

double error(double sum, double sum2, int n) {
  	if (n == 0) return 0;
  	else return sqrt((sum2 - sum*sum) / n);
}

void lattice_RW(int* r, int N_steps, Random &rnd) {
    	for(int k=0; k<N_steps; k++) {
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

void continuum_RW(double* r, int N_steps, Random &rnd) {
	for(int k=0; k<N_steps; k++) {
		double phi = rnd.Rannyu(0, 2*M_PI);
		double theta = rnd.Rannyu(0, M_PI);
		r[0] += sin(theta) * cos(phi);
		r[1] += sin(theta) * sin(phi);
		r[2] += cos(theta); // coordinate sferiche: esercizio mnemonico devastante
	}
}

template<class T> double distance(T* v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

Random init() {
	Random rnd;
   	int seed[4];
   	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
      		Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
      		while ( !input.eof() ){
         		input >> property;
         		if( property == "RANDOMSEED" ){
            			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            			rnd.SetRandom(seed,p1,p2);
         		}
      		}
      		input.close();
   	}
       	else cerr << "PROBLEM: Unable to open seed.in" << endl;
   	return rnd;
}

