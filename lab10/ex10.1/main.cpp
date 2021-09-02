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
#define N_cooling		100
#define N_generations	3000
#define RATE			1.1
#define Starting_BETA	0.05

#include "../../Random/Random.h"
#include "TSP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <numeric>

using namespace std;

vector<City> Set(ofstream&, Random, string);
void TSP(string, ofstream&, ofstream&, ofstream&, Random, double, double);
double Boltzmann(Individual&, double);

int main() {

	// initialize random
	Random rnd;

	// initialize ofstream for both circle and square TSP: the program prints
	// 1 - the cities positions (pos)
	// 2 - the best path between them (best)
	// 3 - the cost of the best individual of the final generation
  	ofstream pos[2];
	ofstream best[2];
	ofstream mean("data/mean.out");
	string type[2] = {"circle", "square"};
	double pmutate[2] = {0.1, 0.25};
	
	// iterate over type = circle and type = square
	for(int i=0; i<2; ++i) {
		pos[i].open("data/" + type[i] + "_pos.out");
		best[i].open("data/" + type[i] + "_best.out");

		// calls TSP algorithm
		TSP(type[i], pos[i], best[i], mean, rnd, Starting_BETA, pmutate[i]);
		pos[i].close();
		best[i].close();
	}

	return 0;
}

/* 
 * This function returns a vector<City> of randomly places
 * cities on a type = "circle" or inside a type = "square".
 */
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

/*
 * This function implements the Traveling Salesman Algorithm
 * with only one individual and then performs metropolis algorithm
 * and finally cools the system every N_generations.
 */
void TSP(string type, ofstream &pos, ofstream& best, ofstream &mean, Random rnd, double beta, double pmutate) {

	// set cities (type = circle or type = square)
	vector<City> cities = Set(pos, rnd, type);

	// initialize vector of ONE individual
	Individual io(cities, rnd);
	vector<Individual> ind;
	ind.push_back(io);

	// start population (only one individual)
	Population p(ind, rnd, N_cities, 0., 0., pmutate);

	for(int i=0; i<N_cooling; ++i) {
		for(int j=0; j<N_generations; ++j) {
			Individual son = p.m_ind[0];
			p.Mutate(son);
			double alpha = Boltzmann(son, beta) / Boltzmann(p.m_ind[0], beta);
			if(alpha >= rnd.Rannyu()) 
				p.m_ind[0] = son;
		}
		mean << p.m_ind[0].Cost() << endl;
		// slowly cool the system (BETA = 1/T, RATE > 1)
		beta *= RATE;
	}

	// print best path
	for(int i=0; i<N_cities; ++i) {
		best << p.m_ind[0].m_city[i].Get_label() << endl;
	}
}

/*
 * This is a typical Boltzmann weight that uses
 * the cost of an individual as its energy
 */
double Boltzmann(Individual &ind, double beta) {
	return exp(-beta * ind.Cost());
}
