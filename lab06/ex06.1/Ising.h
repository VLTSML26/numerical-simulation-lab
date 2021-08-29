/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _Ising_h_
#define _Ising_h_

#include "../../Random/Random.h"

class Ising {

	public:

	Ising();
	Ising(bool, double, double);
	~Ising();

	int Get_step() {return m_nstep;}
	int Get_blk() {return m_nblk;}

	void Reset(int);
	void Accumulate();
	void Averages(int);
	void Move();
	void ConfFinal();
	void Measure();
	double Boltzmann(int, int);
	int Pbc(int);
	void Metropolis();
	void Gibbs();
	void Equilibrate();

	private:

	Random m_rnd;
	int m_props, m_nspin, m_nblk, m_nstep, m_equi;
	bool m_metro;
	double m_beta, m_temp, m_J, m_h;
	double (*m_s), (*m_walker);
	std::ifstream openfile(std::string);

	double blk_av[4], blk_norm;

};

#endif
