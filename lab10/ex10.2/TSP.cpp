/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "TSP.h"
#include <iostream>

using namespace std;

double City::Distance(const City &c) const {
	double x = c.m_x - m_x;
	double y = c.m_y - m_y;

	// returns the L1 distance: hence no sqrt
	return x*x + y*y;
}

double Individual::Cost() const {
	double distance = m_city[m_city.size()].Distance(m_city[0]);
	for(size_t i=0; i<m_city.size()-1; ++i)
		distance += m_city[i].Distance(m_city[i+1]);
	return distance;
}

void Individual::Swap(int a, int b) {
	swap(m_city[a], m_city[b]);
}

bool City::operator==(const City &city) {
	return m_label == city.m_label;
}

bool Individual::operator<(const Individual &ind) {
	return Cost() < ind.Cost();
}

size_t Population::Select() {
	double r = m_rnd.Rannyu(); 
	size_t index = (size_t) (m_ind.size() * pow(r, m_p));
	return index;
}

void Population::Mutate(Individual &ind) {

	// this mutation swaps two random cities
	if(m_rnd.Rannyu() < m_pmutate) {
		int a = (int) (m_rnd.Rannyu(1, m_ncities));
		int b = (int) (m_rnd.Rannyu(1, m_ncities));
		ind.Swap(a, b);
	}

	// this mutation inverts a group of contiguous cities
	if(m_rnd.Rannyu() < m_pmutate) {
		int a = (int) (m_rnd.Rannyu(1, m_ncities));
		int b = (int) (m_rnd.Rannyu(1, m_ncities));
		if(b > a) 
			reverse(ind.m_city.begin() + a, ind.m_city.begin() + b);
		else
			reverse(ind.m_city.begin() + b, ind.m_city.begin() + a);
	}

	// this mutation cycles a group of contiguous cities
	if(m_rnd.Rannyu() < m_pmutate) {
		vector<int> x = {(int) m_rnd.Rannyu(1, m_ncities), (int) m_rnd.Rannyu(1, m_ncities), (int) m_rnd.Rannyu(1, m_ncities)};
		sort(x.begin(), x.end());
		rotate(ind.m_city.begin() + x[0], ind.m_city.begin() + x[1], ind.m_city.begin() + x[2]);
	}
}

void Population::Crossover(Individual p1, Individual p2, Individual &s1, Individual &s2) {
	s1 = p1;
	s2 = p2;
	if(m_rnd.Rannyu() < m_pcross) {
		int cut = (int) (m_rnd.Rannyu(1, m_ncities - 1));
		int index = cut;
		for(int i=0; i<m_ncities; ++i) {
			City c = p2.m_city[i];
			if(find(p1.m_city.begin(), p1.m_city.begin() + cut, c) == p1.m_city.begin() + cut) {
				s1.m_city[index] = c;
				index++;
			}
		}
		index = cut;
		for(int i=0; i<m_ncities; ++i) {
			City c = p1.m_city[i];
			if(find(p2.m_city.begin(), p2.m_city.begin() + cut, c) == p2.m_city.begin() + cut) {
				s2.m_city[index] = c;
				index++;
			}
		}
	}
}

void Population::Evolve() {
	vector<Individual> new_gen;
	for(size_t i=0; i<m_ind.size()/2; ++i) {
		Individual p1 = m_ind[Select()];
		Individual p2 = m_ind[Select()];
		Individual s1, s2;
		Crossover(p1, p2, s1, s2);
		Mutate(s1);
		Mutate(s2);
		new_gen.push_back(s1);
		new_gen.push_back(s2);
	}
	sort(new_gen.begin(), new_gen.end());
	m_ind = new_gen;
}
