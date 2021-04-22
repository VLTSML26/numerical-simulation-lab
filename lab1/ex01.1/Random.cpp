/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Random.h"

using namespace std;

Random::Random() {
  	m_m1 = 502;
  	m_m2 = 1521;
  	m_m3 = 4071;
  	m_m4 = 2107;
  	m_n1 = 0;
  	m_n2 = 0;

	ifstream read = openfile("Primes");
	read >> m_n3 >> m_n4;
	read.close();

	string property;
	read = openfile("seed.in");
	while(!read.eof()) {
		read >> property;
		if(property == "RANDOMSEED")
			read >> m_l1 >> m_l2 >> m_l3 >> m_l4;
	}
	read.close();

  	m_l1 = m_l1 % 4096;
  	m_l2 = m_l2 % 4096;
  	m_l3 = m_l3 % 4096;
  	m_l4 = m_l4 % 4096;
  	m_l4 = 2 * (m_l4/2) + 1;
}

void Random::SaveSeed() {
   	ofstream WriteSeed;
   	WriteSeed.open("seed.out");
   	if (WriteSeed.is_open()) {
      	WriteSeed << m_l1 << "\t" << m_l2 << "\t" << m_l3 << "\t" << m_l4 << endl;;
	} 
	else 
		cerr << "PROBLEM: Unable to open random.out" << endl;
  	WriteSeed.close();
  	return;
}

double Random::Gauss(double mean, double sigma) {
   	double s = Rannyu();
   	double t = Rannyu();
   	double x = sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   	return mean + x * sigma;
}

double Random::Exp(double lambda) {
	double s = Rannyu();
	double x = - log(1. - s);
	return x / lambda;
}

double Random::Lorentz(double mu, double gamma) {
	double s = Rannyu();
	s -= 0.5;
	s *= M_PI;
	double x = tan(s);
	return mu + x * gamma;
}

double Random::Line() {
	double s = Rannyu();
	return 1 - sqrt(1 - s);
}

bool Random::Bool() {
	double s = Rannyu(0,2);

	if((int)(s) == 0) 
		return true;
	else 
		return false;
}

double Random::Rannyu(double min, double max) {
   	return min + (max - min) * Rannyu();
}

double Random::Rannyu(void) {
  	const double twom12 = 0.000244140625;

  	int i1 = m_l1*m_m4 + m_l2*m_m3 + m_l3*m_m2 + m_l4*m_m1 + m_n1;
  	int i2 = m_l2*m_m4 + m_l3*m_m3 + m_l4*m_m2 + m_n2;
  	int i3 = m_l3*m_m4 + m_l4*m_m3 + m_n3;
  	int i4 = m_l4*m_m4 + m_n4;

  	m_l4 = i4 % 4096;
  	i3 = i3 + i4/4096;
  	m_l3 = i3 % 4096;
  	i2 = i2 + i3/4096;
  	m_l2 = i2 % 4096;
  	m_l1 = (i1 + i2/4096) % 4096;

  	double r = twom12 * (m_l1 + twom12*(m_l2 + twom12*(m_l3 + twom12*(m_l4))));

  	return r;
}

ifstream Random::openfile(string file) {
	ifstream read;
	read.open(file);
	if(!read.is_open()) {
		cerr << "Error: unable to open " << file << endl;
		exit(0);
	}
	return read;
}
