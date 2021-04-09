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

  	const int N_throws = 10000;
  	const int N_blocks = 100;

  	ofstream data1 ("ex02.1.1.out"); // uniform sampling
  	ofstream data2 ("ex02.1.2.out"); // p(x) = 2(1-x)

	// ==== EX02.1.1: UNIFORM SAMPLING ====
	
	double sum_prog = 0., sum2_prog = 0.;
   	
	for(int i=0; i<N_blocks; i++){
     		double sum = 0.;
     		for(int j=0; j<N_throws; j++){
       			double ran = rnd.Rannyu();
       			double y = M_PI / 2 * cos(M_PI / 2 * ran);
			sum += y;
     		}
		// ==== ANALISI DATI (metodo a blocchi) ====
     		sum_prog += sum / N_throws;
     		sum2_prog += sum / N_throws * sum / N_throws;
     		double av = sum_prog / (i + 1);
     		double av2 = sum2_prog / (i + 1);
     		data1 << av << "\t" << error(av, av2, i) << endl;
	}

	// ==== EX02.1.2: SAMPLING CON p(x) = 2(1-x) ====
	
	sum_prog = 0.; sum2_prog = 0.;
   
	for(int i=0; i<N_blocks; i++){
     		double sum = 0.;
     		for(int j=0; j<N_throws; j++){
       			double ran = rnd.Ex02_1();
       			double y = (M_PI / 2 * cos(M_PI / 2 * ran)) / (2 * (1 - ran));
			sum += y;
     		}
		// ==== ANALISI DATI (metodo a blocchi) ====
     		sum_prog += sum / N_throws;
     		sum2_prog += sum / N_throws * sum / N_throws;
     		double av = sum_prog / (i + 1);
     		double av2 = sum2_prog / (i + 1);
     		data2 << av << "\t" << error(av, av2, i) << endl;
	}
	
	rnd.SaveSeed();
   	return 0;
}

double error(double sum, double sum2, int n) {
  	if (n == 0) return 0;
  	else return sqrt((sum2 - sum*sum) / n);
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
