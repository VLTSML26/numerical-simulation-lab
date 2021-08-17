/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _TSP_h_
#define _TSP_h_

#include "../../Random/Random.h"
#include <string>
#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>

class City {

	public :

	City() {m_x = 0; m_y = 0;}
	City(double x, double y) {m_x = x; m_y = y;}
	~City() {;}

	double Distance(City);

	double Get_x() {return m_x;}
	double Get_y() {return m_y;}
	int Get_label() {return m_label;}

	private :

	double m_x, m_y;
	int m_label;

};

class Individual {

	public :
	
	Individual(std::vector<City> city) {m_city = city;}

	double Cost();
	void Swap(int, int);
	bool operator<(Individual&);

	private :

	std::vector<City> m_city;

};

class Population {

	public :

	Population(std::vector<Individual> ind) {m_ind = ind;}

	void Evolve();

	private :

	std::vector<Individual> m_ind;
	Individual Mutate(Individual);
	int Select();
	int m_ncities;
	Random m_rnd;
	double m_p, m_pcross, m_pmutate;

};

#endif
