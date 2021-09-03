/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Ising.h"

using namespace std;

/*
 * Constructor
 */
Ising::Ising(bool metropolis, double t, double h) {

	m_props = 4;
	m_metro = metropolis;
	m_temp = t;
	m_h = h;

	// reads config data from input.dat
	ifstream read = openfile("input/input.dat");
	read >> m_nspin >> m_J >> m_nblk >> m_nstep >> m_equi;
	m_beta = 1. / m_temp;
	read.close();

	m_s = new double[m_nspin];
	ifstream config;
	config.open("input/config.final");
	if(!config.is_open()) {
		for(int i=0; i<m_nspin; ++i) {
			if(m_rnd.Rannyu() >= 0.5) 
				m_s[i] = 1;
			else 
				m_s[i] = -1;
		}
	} else {
		m_rnd = Random("../../Random/seed.out");
		for(int i=0; i<m_nspin; ++i)
			config >> m_s[i];
	}

	// 0 - Energy
	// 1 - Heat Capacity
	// 2 - Magnetic susceptibility
	// 3 - Magnetization
	m_walker = new double[m_props];
	for(int i=0; i<m_props; ++i)
			m_walker[i] = 0.;
}

Ising::~Ising() {
	delete m_s;
	delete m_walker;
}

/*
 * This method uses Metropolis algorithm to sample the Boltzmann weight
 */
void Ising::Metropolis() {
	int i = (int) m_rnd.Rannyu(0, m_nspin);
	double alpha = exp(m_beta * Boltzmann(m_s[i], i));
	if(alpha >= m_rnd.Rannyu())
		m_s[i] = - m_s[i];
}

/*
 * This method uses Gibbs algorithm to sample the Boltzmann weight
 */
void Ising::Gibbs() {
	int i = (int) m_rnd.Rannyu(0, m_nspin);
	m_s[i] = 1;
	double alpha = exp(m_beta * Boltzmann(m_s[i], i));
	if(m_rnd.Rannyu() > 1. / (1. + alpha))
		m_s[i] = -1;
}

void Ising::Move() {
	for(int i=0; i<m_nspin; ++i) {
		if(m_metro)
			Metropolis();
		else
			Gibbs();
	}
}

/*
 * This is the Boltzmann weight
 */
double Ising::Boltzmann(int sm, int ip) {
	double ene = - m_J * sm * (m_s[Pbc(ip - 1)] + m_s[Pbc(ip + 1)]) - m_h * sm;
	return 2 * ene;
}

void Ising::Measure() {
	double u = 0., s = 0.;

	for (int i=0; i<m_nspin; ++i) {
		u += - m_J * m_s[i] * m_s[Pbc(i+1)] - 0.5 * m_h * (m_s[i] + m_s[Pbc(i+1)]);
		s += m_s[i];
	}

	m_walker[0] = u;
	m_walker[1] = u * u;
	m_walker[2] = m_beta * s * s;
	m_walker[3] = s;
}

/*
 * Data blocking 
 */
void Ising::Reset(int iblk) {
   
	for(int i=0; i<m_props; ++i)
		blk_av[i] = 0;

	blk_norm = 0;
}

/*
 * Data blocking 
 */
void Ising::Accumulate() {

	for(int i=0; i<m_props; ++i)
		blk_av[i] = blk_av[i] + m_walker[i];

	blk_norm = blk_norm + 1.0;
}

/*
 * Data blocking 
 */
void Ising::Averages(int iblk) {
    
	string name;
	if(m_h == 0) {
		if(m_metro)
			name = "metropolis.out";
		else
			name = "gibbs.out";
	} else {
		if(m_metro)
			name = "metropolis_extfield.out";
		else
			name = "gibbs_extfield.out";
	}
	ofstream write;
	write.open("data/" + name, ios::app);
	
	double e = blk_av[0] / blk_norm / (double) m_nspin;
	double c = m_beta * m_beta * (blk_av[1] / blk_norm / (double) m_nspin - (double) m_nspin * e * e);
	double chi = blk_av[2] / blk_norm / (double) m_nspin;
	double M = blk_av[3] / blk_norm / (double) m_nspin;

	write << e << "\t" << c << "\t" << chi << "\t" << M << endl;
    write.close();
}

/*
 * This method prints out the final spin configuration
 */
void Ising::ConfFinal() {

	ofstream write;
	write.open("input/config.final");

	for (int i=0; i<m_nspin; ++i)
		write << m_s[i] << endl;

	write.close();
	m_rnd.SaveSeed();
}

/*
 * This method implements periodic boundary conditions
 */
int Ising::Pbc(int i) {
    if(i >= m_nspin) 
		i -= m_nspin;
    else if(i < 0) 
		i += m_nspin;
    return i;
}

/*
 * This method usually called at the beginning of the simulation
 * equilibrates the system (ie the Metropolis/Gibbs algorithm)
 */
void Ising::Equilibrate() {
	for(int i=0; i<m_equi; ++i)
		Move();
}

ifstream Ising::openfile(string file) {
	ifstream read;
	read.open(file);
	if(!read.is_open()) {
		cerr << "Error: unable to open " << file << endl;
		exit(0);
	}
	return read;
}
