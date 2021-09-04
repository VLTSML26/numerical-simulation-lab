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

/*
 * This class represents a city of coordinates (x,y)
 * and integer label.
 */
class City {

	public :

	City() {m_x = 0; m_y = 0; m_label = 0;}
	City(double x, double y, int label) {m_x = x; m_y = y; m_label = label;}
	~City() {;}

	// returns the L1 distance between two cities
	double Distance(const City&) const;

	double Get_x() {return m_x;}
	double Get_y() {return m_y;}
	int Get_label() {return m_label;}

	// overloading of operator==, useful with std::vector
	bool operator==(const City&);

	private :

	double m_x, m_y;
	int m_label;

};

/*
 * This class represents an individual as 
 * and ordered path of cities
 */
class Individual {

	public :
	
	Individual() {;}

	// constructor: randomly stores the cities from input
	// only maintaining the first one in that position 
	Individual(std::vector<City> c, Random &rnd) {
		m_city.push_back(c[0]);
		for(size_t i=1; i<c.size(); ++i) {
			for(;;) {
				City appo = c[rnd.Rannyu(1., c.size())];
				if(std::find(m_city.begin() + 1, m_city.end(), appo) == m_city.end()) {
					m_city.push_back(appo);
					break;
				}
			}
		}
	} 

	// returns the cost: i.e. the distance covered by visiting
	// the cities in the order set for this individual
	double Cost() const;

	// swaps cities
	void Swap(int, int);

	// overloading of operator<, useful with std::vector
	bool operator<(const Individual&);

	// this member used to be private, but it is required outside the class
	// a huge number of times
	std::vector<City> m_city;

};

/*
 * This class represents a population as a std::vector of individuals
 * and contains its evolution algorithms
 */
class Population {

	public :

	Population(std::vector<Individual> ind, Random &rnd, int cities, double p, double pcross, double pmutate) :	
		m_ind{ind}, m_rnd{rnd},	m_ncities{cities}, m_p{p}, m_pcross{pcross}, m_pmutate{pmutate} 
		{
			std::sort(m_ind.begin(), m_ind.end());
		}

	// evolves the population
	void Evolve();

	// crossover happens with probability m_pcross
	void Crossover(Individual, Individual, Individual&, Individual&);

	// mutates an individual with probability m_pmutate
	void Mutate(Individual&);

	// sort the population
	void Sort() {
		std::sort(m_ind.begin(), m_ind.end());
	}

	std::vector<Individual> m_ind;

	private :

	Random &m_rnd;
	size_t Select();
	int m_ncities;
	double m_p, m_pcross, m_pmutate;

};

#endif
