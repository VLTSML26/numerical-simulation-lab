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
#include "TSP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

int main() {

	Random rnd;

	int N_cities = 32;
	int N_individuals = 5000;
	int N_generations = 300;
	
	vector<City> cities;
	vector<Individual> ind;

	for(int i=0; i<N_cities; ++i) {
		double theta = rnd.Rannyu(0, 2*M_PI);
		City c(cos(theta), sin(theta), i);
		cities.push_back(c);
	}

	for(int i=0; i<N_individuals; ++i) {
		Individual io(cities, rnd);
		ind.push_back(io);
	}

	Population p(ind, rnd);

	for(int i=0; i<N_generations; ++i)
		p.Evolve();

	return 0;
}
