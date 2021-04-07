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
  	int N_path = 100; // we divide the interval [0, t1] in N_path equal bins and we sample the discretized GMB path of the asset price
	double s0 = 100.; // asset price at t=0
	double k = 100.; // strike price
	double sigma = 0.25; // volatility
	double r = 0.1; // risk-free interest rate
	double t1 = 1.; // delivery time

  	ofstream data1 ("ex03.1.1.out");
	ofstream data2 ("ex03.1.2.out");

	// ==== NOTA BENE ====
	// Le variabili che terminano con "_prog" vengono inizializzate qui sotto e non devono MAI essere toccate all'interno dei cicli

	// ==== EX03.1: DIRECT SAMPLE ====
	// C stands for "Call Option"
	// P stands for "Put Option"

	double C_sum_prog = 0., C_sum2_prog = 0.;
	double P_sum_prog = 0., P_sum2_prog = 0.;
   	for(int i=0; i<N_blocks; i++){
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; j++){
			double z = rnd.Gauss(0., 1.);
			double s = s0 * exp((r - sigma*sigma/2.)*t1 + sigma*z*sqrt(t1));
			s -= k; // NOTA BENE
			if(s > 0) C_sum += s*exp(- r*t1); // il profitto per una Call è max(0, s - k) ovvero max(0, s) vedi riga sopra
			else P_sum += - s*exp(- r*t1); // il profitto per una Put è max(0, k - s) ovvero max(0, -s) ibidem
		}
		// Call data analysis
		C_sum_prog += C_sum / N_throws;
     	C_sum2_prog += C_sum / N_throws * C_sum / N_throws;
     	double C_av = C_sum_prog / (i + 1);
     	double C_av2 = C_sum2_prog / (i + 1);
     	// Put data analysis
     	P_sum_prog += P_sum / N_throws;
     	P_sum2_prog += P_sum / N_throws * P_sum / N_throws;
     	double P_av = P_sum_prog / (i + 1);
     	double P_av2 = P_sum2_prog / (i + 1);
     	data1 << C_av << "\t" << error(C_av, C_av2, i) << "\t" << P_av << "\t" << error(P_av, P_av2, i) << endl;
	}

	// ==== EX03.2: DISCRETIZED SAMPLE ====
	// C stands for "Call Option"
	// P stands for "Put Option"

	C_sum_prog = 0., C_sum2_prog = 0.;
	P_sum_prog = 0., P_sum2_prog = 0.;
	
   	for(int i=0; i<N_blocks; i++){
     	double C_sum = 0.;
     	double P_sum = 0.;
     	for(int j=0; j<N_throws; j++){
			double s = s0;
			for(int k=0; k<N_path; k++) {
				double z = rnd.Gauss(0., 1.);
				s *= exp((r - sigma*sigma/2.)*t1/N_path + sigma*z*sqrt(t1/N_path));
			}
			s -= k; // NOTA BENE
			if(s > 0) C_sum += s*exp(- r*t1); // il profitto per una Call è max(0, s - k) ovvero max(0, s) vedi riga sopra
			else P_sum += - s*exp(- r*t1); // il profitto per una Put è max(0, k - s) ovvero max(0, -s) ibidem

		}
		// Call data analysis
		C_sum_prog += C_sum / N_throws;
     	C_sum2_prog += C_sum / N_throws * C_sum / N_throws;
     	double C_av = C_sum_prog / (i + 1);
     	double C_av2 = C_sum2_prog / (i + 1);
     	// Put data analysis
     	P_sum_prog += P_sum / N_throws;
     	P_sum2_prog += P_sum / N_throws * P_sum / N_throws;
     	double P_av = P_sum_prog / (i + 1);
     	double P_av2 = P_sum2_prog / (i + 1);
     	data2 << C_av << "\t" << error(C_av, C_av2, i) << "\t" << P_av << "\t" << error(P_av, P_av2, i) << endl;
	}

	rnd.SaveSeed();
   	return 0;
}

double error(double D_sum, double D_sum2, int n) {
  	if (n == 0) return 0;
  	else return sqrt((D_sum2 - D_sum*D_sum) / n);
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