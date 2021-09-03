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
#include "NVT.h"

using namespace std;

NVT::NVT(string folder) : MolDyn(folder) {

	m_nbins = 100;
	m_accepted = 0;
	m_props = 2 + m_nbins;
	m_beta = 1. / m_temp;
	m_binsize = (m_box / 2.) / (double) m_nbins;
	
	// tail corrections
	m_vtail = (8. * M_PI * m_rho) / (9. * pow(m_rcut, 9)) - (8. * M_PI * m_rho) / (3. * pow(m_rcut, 3));
	m_ptail = (32. * M_PI * m_rho) / (9. * pow(m_rcut, 9)) - (16. * M_PI * m_rho) / (3. * pow(m_rcut, 3));
	
	// 0 - Potential Energy
	// 1 - Virial
	// 2 - 1st bin g(r)
	// ...
	// 2 + m_nbin - last bin g(r)
	m_walker = new double[m_props];
	for(int i=0; i<m_props; ++i)
			m_walker[i] = 0.;
}

NVT::~NVT() {
	delete m_walker;
}

void NVT::Move() {
	double xold, yold, zold, xnew, ynew, znew;
	
	for(int i=0; i<m_npart; ++i) {
		int k = (int) (m_rnd.Rannyu() * m_npart);
	
		xold = m_x[k];
		yold = m_y[k];
		zold = m_z[k];
	
		double energy_old = Boltzmann(xold, yold, zold, k);

		xnew = Pbc(m_x[k] + m_dt * (m_rnd.Rannyu() - 0.5));
		ynew = Pbc(m_y[k] + m_dt * (m_rnd.Rannyu() - 0.5));
		znew = Pbc(m_z[k] + m_dt * (m_rnd.Rannyu() - 0.5));
	
		double energy_new = Boltzmann(xnew,ynew,znew,k);
	
		double p = exp(m_beta*(energy_old - energy_new));
		if(p >= m_rnd.Rannyu()) {
			m_x[k] = xnew;
			m_y[k] = ynew;
			m_z[k] = znew;
			m_accepted++;
		}
	}
}

void NVT::Tune(int n) {
	bool half = true;
	double rate = 0;
	do {
		m_accepted = 0;
		for(int i=0; i<n; ++i)
			Move();
		rate = (double) m_accepted / (n * m_npart);
		cout << m_dt << "\t" << rate << endl;
		if(rate > 0.62)
			m_dt *= 1.1;
		else if(rate < 0.38)
			m_dt *= 0.9;
		else half = false;
	} while(half);
	m_accepted = 0;
	cout << "Tuned with delta: " << m_dt << ". Acceptance rate: " << rate << endl;
}

double NVT::Boltzmann(double xx, double yy, double zz, int ip) {
	double ene = 0.;
	for (int i=0; i<m_npart; ++i) {
		if(i != ip) {
			double dx = Pbc(xx - m_x[i]);
			double dy = Pbc(yy - m_y[i]);
			double dz = Pbc(zz - m_z[i]);
		
			double dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
		
			if(dr < m_rcut) {
				ene += 1. / pow(dr,12) - 1. / pow(dr,6);
			}
		}
	}
	return 4.0 * ene;
}

void NVT::Measure() {
	double v = 0., w = 0.;
	
	for(int k=0; k<m_props; ++k)
		m_walker[k] = 0.;
	
	for(int i=0; i<m_npart-1; ++i) {
		for(int j=i+1; j<m_npart; ++j) {
			double dx = Pbc(m_x[i] - m_x[j]);
			double dy = Pbc(m_y[i] - m_y[j]);
			double dz = Pbc(m_z[i] - m_z[j]);
			double dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			// g(r)
			if(dr < m_box / 2.) {
				int bin = (int) (2 * dr * m_nbins / m_box);
				m_walker[2 + bin] += 2;
			}
			// properties
			if(dr < m_rcut) {
				double vij = 1. / pow(dr, 12) - 1. / pow(dr, 6);
				double wij = 1. / pow(dr, 12) - 0.5 / pow(dr, 6);
				v += vij;
				w += wij;
			}
		}          
	}

	m_walker[0] = 4. * v;
	m_walker[1] = 48. * w / 3.;
}

/*
 * This method is usually called after Measure() in order to
 * add tail corrections and works also as a rescale of the bins'
 * occupation number for a correct estimation of g(r)
 */
void NVT::Measure(double sum[]) {
	double vol = (double) m_npart / m_rho;
	sum[0] = sum[0] / m_throws / (double) m_npart + m_vtail;
    sum[1] = m_rho * m_temp + (sum[1] / m_throws + m_ptail * (double) m_npart) / vol;
	for(int i=0; i<m_nbins; i++) {
		double bin_size = (m_box / 2.0) / (double)m_nbins;
		double r = i * bin_size;
		double deltaV = 4. * M_PI/3. * (pow((r+bin_size), 3.) - pow(r, 3.));
		double div = m_rho * m_npart * deltaV * m_throws;
		sum[2 + i] /= div;
	}
}

void NVT::ConfFinal() {
	ofstream write;
	write.open("config.final");
	for(int i=0; i<m_npart; ++i)
		write << m_x[i] / m_box << "\t" << m_y[i] / m_box << "\t" << m_z[i] / m_box << endl;
	write.close();
	
	m_rnd.SaveSeed();
}
