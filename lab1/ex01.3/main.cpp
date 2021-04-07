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
Random init();

int main (int argc, char *argv[]){

	Random rnd;
	rnd = init();

  	int N_throws = 10000;
  	int N_blocks = 100;
  	double d = 1., l = 0.7;

  	ofstream data ("ex01.3.1.out");

	// ==== NOTA BENE ====
	// Le variabili che terminano con "_prog" vengono inizializzate qui sotto e non devono MAI essere toccate all'interno dei cicli	

	double sum_prog = 0., sum2_prog = 0.;
   	for(int i=0; i<N_blocks; i++){
     	int hit = 0;
     	for(int j=0; j<N_throws; j++){
       		double x1 = rnd.Rannyu(0., d);
       		double x, y, theta;
       		do {
       			x = rnd.Rannyu();
       			y = rnd.Rannyu();
       		} while (x*x + y*y > 1.);
       		if(y >= 0.) theta = acos(x/sqrt(x*x + y*y));
       		if(y < 0.) theta = 2*M_PI - acos(x/sqrt(x*x + y*y));
       		double x2 = x1 + l*sin(theta);
       		if(x2 < 0. or x2 > d) hit++;
     	}
	// ==== ANALISI DATI (metodo a blocchi) ====
		double pi = 2*l*N_throws / (hit*d);
     	sum_prog += pi;
     	sum2_prog += pi * pi;
     	double av = sum_prog / (i + 1);
     	double av2 = sum2_prog / (i + 1);
     	data << av << "\t" << error(av, av2, i) << endl;
	}	

	rnd.SaveSeed();
   	return 0;
}

double error(double sum, double sum2, int n) {
  	if (n == 0) {
    	return 0;
    }
  	else  {
    	return sqrt((sum2 - sum*sum) / n);
    }
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