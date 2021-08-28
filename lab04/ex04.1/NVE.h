/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _NVE_h_
#define _NVE_h_

#include "../../MolDyn/MolDyn.h"

class NVE : public MolDyn {

	public :

	NVE(bool, std::string old = std::string());
   	~NVE();

	virtual void Move();
	virtual void Measure();

	double Force(int, int);
	void ConfFinal(std::string, std::string);
	double walker(int k) {return m_walker[k];}

	private :

	double (*m_xold), (*m_yold), (*m_zold), // xold
		   (*m_vx), (*m_vy), (*m_vz); // v

};

#endif
