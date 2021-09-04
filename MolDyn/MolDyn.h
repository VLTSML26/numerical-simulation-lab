/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _Mol_Dyn_h_
#define _Mol_Dyn_h_

#include "../Random/Random.h"
#include <string>
#include <fstream>

class MolDyn {
	
	public :
	
	MolDyn(std::string);
	~MolDyn();

	// virtual methods will be implemented differently
	// by NVE or NVT
	virtual void Move() = 0;
	virtual void Measure() = 0;
	
	int Get_throws() {return m_throws;}
	int Get_npart() {return m_npart;}
	int Get_props() {return m_props;}
	int Get_blocks() {return m_blocks;}
	double Pbc(double);

	// walker is the array that contains measured observables
	double walker(int k) {return m_walker[k];}
	
	// this is a virtual class, no private members
	protected :

	Random m_rnd;
	int m_npart, m_blocks, m_throws, m_props, m_nbins;
	double m_temp, m_rcut, m_dt, m_box, m_rho;
	double (*m_x), (*m_y), (*m_z), (*m_walker);

	std::ifstream openfile(std::string);

};

#endif
