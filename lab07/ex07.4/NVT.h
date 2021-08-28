/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _NVT_h_
#define _NVT_h_

#include "../../Random/Random.h"
#include "../../MolDyn/MolDyn.h"

class NVT : public MolDyn {

	public :

	NVT();
	~NVT();

	virtual void Move();
	virtual void Measure();

	void ConfFinal();
	void Measure(double[]);
	double Boltzmann(double, double, double, int);

	private :

	double m_beta, m_vtail, m_ptail, m_binsize;
};

#endif
