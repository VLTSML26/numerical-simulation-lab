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
#define N_individuals	5000
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
void TSP(string, ofstream&, ofstream&, ofstream&, Random, double, double, double);
double Accumulate(vector<Individual>);

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

	// set the probabilities for each process
	double exp[2] = {3., 4.};
	double pcross[2] = {0.75, 0.6};
	double pmutate[2] = {0.55, 0.05};
	
	// iterate over type = circle and type = square
	for(int i=0; i<2; ++i) {
		pos[i].open("data/" + type[i] + "_pos.out");
		best[i].open("data/" + type[i] + "_best.out");

		// calls TSP algorithm
		TSP(type[i], pos[i], best[i], mean, rnd, exp[i], pcross[i], pmutate[i]);
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
 * and cools the system every N_generations.
 */
void TSP(string type, ofstream &pos, ofstream& best, ofstream &mean, Random rnd, double exp, double pcross, double pmutate) {

	// set cities (type = circle or type = square)
	vector<City> cities = Set(pos, rnd, type);

	// initialize vector of N_individuals individuals
	vector<Individual> ind;
	for(int i=0; i<N_individuals; ++i) {							
		Individual io(cities, rnd);
		ind.push_back(io);
	}

	// generate population
	Population p(ind, rnd, N_cities, exp, pcross, pmutate);

	// evolve population for N_generations times
	for(int i=0; i<N_generations; ++i) {
		p.Evolve();
		mean << Accumulate(p.m_ind) << "\t" << p.m_ind[0].Cost() << endl;
	}

	for(int i=0; i<N_cities; ++i) {
		best << p.m_ind[0].m_city[i].Get_label() << endl;
	}
}

double Accumulate(vector<Individual> ind) {
	return accumulate(ind.begin(), ind.begin() + N_individuals / 2, ind[0].Cost(), [](const double x, const Individual &io) {return x + io.Cost();}) / (N_individuals / 2);
}
