/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "MolDyn.h"
#include <cmath>
#include <iostream>

using namespace std;

MolDyn::MolDyn() {

	ifstream read = openfile("../input/input.dat");
	read >> m_temp >> m_npart >> m_rho >> m_rcut >> m_dt >> m_blocks >> m_throws;
	read.close();

	double vol = (double) m_npart / m_rho;
	m_box = pow(vol, 1./3.);

	m_x = new double[m_npart];
	m_y = new double[m_npart];
	m_z = new double[m_npart];

	read = openfile("../input/config.0");
  	for(int i=0; i<m_npart; ++i) {
    	read >> m_x[i] >> m_y[i] >> m_z[i];
    	m_x[i] *= m_box;
    	m_y[i] *= m_box;
    	m_z[i] *= m_box;
  	}
  	read.close();
}

MolDyn::~MolDyn() {
	m_rnd.SaveSeed();
	delete m_x;
	delete m_y;
	delete m_z;
	delete m_walker;
}

double MolDyn::Pbc(double r) { 
    	return r - m_box * rint(r / m_box);
}

ifstream MolDyn::openfile(string file) {
	ifstream read;
	read.open(file);
	if(!read.is_open()) {
		cerr << "Error: unable to open " << file << endl;
		exit(0);
	}
	return read;
}
