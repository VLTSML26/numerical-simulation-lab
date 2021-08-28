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

	City() {m_x = 0; m_y = 0; m_label = 0;}
	City(double x, double y, int label) {m_x = x; m_y = y; m_label = label;}
	~City() {;}

	double Distance(const City&) const;

	double Get_x() {return m_x;}
	double Get_y() {return m_y;}
	int Get_label() {return m_label;}
	bool operator==(const City&);

	private :

	double m_x, m_y;
	int m_label;

};

class Individual {

	public :
	
	Individual() {;}
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

	double Cost() const;
	void Swap(int, int);
	bool operator<(const Individual&);

	std::vector<City> m_city;

};

class Population {

	public :

	Population(std::vector<Individual> ind, Random &rnd, int cities) :	m_ind{ind}, m_rnd{rnd},	m_ncities{cities} 
																		{
																			m_p = 4.; 
																			m_pcross = 0.75;
																			m_pmutate = 0.5;
																			std::sort(m_ind.begin(), m_ind.end());
																		}

	void Evolve();
	void Crossover(Individual, Individual, Individual&, Individual&);
	void Mutate(Individual&);
	std::vector<Individual> m_ind;

	private :

	Random &m_rnd;
	size_t Select();
	int m_ncities;
	double m_p, m_pcross, m_pmutate;

};

#endif
