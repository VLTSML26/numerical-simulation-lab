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
 
Random init();

int main (int argc, char *argv[]){

	Random rnd;
	rnd = init();

  	int N_throws = 10000;
  	const int Ns[4] = {1, 2, 10, 100};
  	const char stampa[4] = {'\t', '\t', '\t', '\n'};

  	ofstream data1 ("ex01.2.1.out"); // uniform
  	ofstream data2 ("ex01.2.2.out"); // exponential
  	ofstream data3 ("ex01.2.3.out"); // lorentz

	// ==== NOTA BENE ====
	// Le variabili che terminano con "_prog" vengono inizializzate qui sotto e non devono MAI essere toccate all'interno dei cicli	
	// U sta per Uniform
	// E sta per Exponential
	// L sta per Lorentz
	
   	for(int i=0; i<N_throws; i++){
     	for(int j=0; j<4; j++){
       		double U_sum = 0;
       		double E_sum = 0;
       		double L_sum = 0;
       		for(int k=0; k<Ns[j]; k++) {
       			U_sum += (int) rnd.Rannyu(1., 7.);
       			E_sum += rnd.Exp(1.);
       			L_sum += rnd.Lorentz(0., 1.);
       		}
       		data1 << (double) (U_sum / Ns[j]) << stampa[j];
       		data2 << E_sum / Ns[j] << stampa[j];
       		data3 << L_sum / Ns[j] << stampa[j];
     	}
	}	

	rnd.SaveSeed();
   	return 0;
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
   	string property; //D=26
   	if (input.is_open()){
      	while ( !input.eof() ){
         	input >> property;
         	if( property == "RANDOMSEED" ){
            	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            	rnd.SetRandom(seed,p1,p2);
         	}
      	}
      	input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
   	return rnd;
}