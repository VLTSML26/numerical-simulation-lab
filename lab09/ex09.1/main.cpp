/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_cities		32
#define N_individuals	6000
#define N_generations	350

#include "../../Random/Random.h"
#include "TSP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <numeric>

using namespace std;

vector<City> Set(ofstream&, Random, string);
void TSP(string, ofstream&, ofstream&, ofstream&, Random);
double Accumulate(vector<Individual>);

int main() {

	Random rnd;													// initialize Random
  	ofstream pos[2];											// output stream 
	ofstream best[2];
	ofstream mean("mean.out");
	string type[2] = {"circle", "square"};							// 
	
	for(int i=0; i<2; ++i) {
		pos[i].open(type[i] + "_pos.out");
		best[i].open(type[i] + "_best.out");
		TSP(type[i], pos[i], best[i], mean, rnd);
		pos[i].close();
		best[i].close();
	}

	return 0;
}

vector<City> Set(ofstream &pos, Random rnd, string type) {
	vector<City> cities;

	if(type == "circle") {
		for(int i=0; i<N_cities; ++i) {
			double theta = rnd.Rannyu(0, 2*M_PI);
			City c(cos(theta), sin(theta), i);
			cities.push_back(c);
			pos	<< c.Get_x() << "\t" << c.Get_y() << endl;
		}
	}
	else if(type == "square") {
		for(int i=0; i<N_cities; ++i) {
			double x = rnd.Rannyu(0, 1);
			double y = rnd.Rannyu(0, 1);
			City c(x, y, i);
			cities.push_back(c);
			pos	<< c.Get_x() << "\t" << c.Get_y() << endl;
		}
	}
	else {
		cout << "Error: for now, let's place the cities on a <circle> or a <square>" << endl;
		exit(0);
	}

	return cities;
}

void TSP(string type, ofstream &pos, ofstream& best, ofstream &mean, Random rnd) {
	vector<City> cities = Set(pos, rnd, type);						// see function above

	vector<Individual> ind;											// initialize vector of individuals to start a population
	for(int i=0; i<N_individuals; ++i) {							
		Individual io(cities, rnd);									// fill vector of individuals
		ind.push_back(io);
	}

	Population p(ind, rnd, N_cities);								// start population

	for(int i=0; i<N_generations; ++i) {
		p.Evolve();													// evolve population N_generatons times
		mean << Accumulate(p.m_ind) << "\t" << p.m_ind[0].Cost() << endl;
	}

	for(int i=0; i<N_cities; ++i) {
		best << p.m_ind[0].m_city[i].Get_label() << endl;
	}
}

double Accumulate(vector<Individual> ind) {
	return accumulate(ind.begin(), ind.begin() + N_individuals / 2, ind[0].Cost(), [](const double x, const Individual &io) {return x + io.Cost();}) / (N_individuals / 2);
}
