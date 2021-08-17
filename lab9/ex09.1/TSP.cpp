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
}

void Population::Crossover(Individual p1, Individual p2, Individual &s1, Individual &s2) {
	int cut = (int) (1 + (m_ncities - 1) * m_rnd.Rannyu());
	int index = cut;
	for(int i=0; i<m_ncities; ++i) {
		if(find(

void Population::Evolve() {
	vector<Individual> new_gen;
	for(int i=0; i<m_ind.size(); ++i) {
		// p stands for parent
		Individual p1 = m_ind[Select()];
		Individual p2 = m_ind[Select()];
	//	if(m_rnd.Rannyu() < m_pcross) Crossover();
		Individual son1 = Mutate(p1);
		Individual son2 = Mutate(p2);
	}
	//sort(m_ind)
}

