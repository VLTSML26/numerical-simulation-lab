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

using namespace std;

double City::Distance(City c) {
	double x = c.m_x - m_x;
	double y = c.m_y - m_y;
	return x*x + y*y; // tolto sqrt
}

double Individual::Cost() {
	double distance = 0.;
	for(auto &c : m_city)
		distance += c.Distance(*(&c+1));
	return distance;
}

void Individual::Swap(int a, int b) {
	swap(m_city[a], m_city[b]);
}

bool City::operator==(const City &city) {
	return m_label == city.m_label;
}

bool Individual::operator<(Individual &ind) {
	return Cost() < ind.Cost();
}

int Population::Select() {
	double r = m_rnd.Rannyu(); 
	int index = (int) (m_ind.size() * pow(r, m_p));
	return index;
}

Individual Population::Mutate(Individual ind) {
	if(m_rnd.Rannyu() < m_pmutate) {
		int a = (int) (1 + (m_ncities - 1) * m_rnd.Rannyu());
		int b = (int) (1 + (m_ncities - 1) * m_rnd.Rannyu());
		ind.Swap(a, b);
	}
	return ind;
}

void Population::Crossover(Individual p1, Individual p2, Individual &s1, Individual &s2) {
	s1 = p1;
	s2 = p2;
	if(m_rnd.Rannyu() < m_pcross) {
		int cut = (int) (1 + (m_ncities - 1) * m_rnd.Rannyu());
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
		s1 = Mutate(s1);
		s2 = Mutate(s2);
		new_gen.push_back(s1);
		new_gen.push_back(s2);
	}
	sort(new_gen.begin(), new_gen.end());
	m_ind = new_gen;
}
