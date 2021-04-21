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

#include <string>
#include <fstream>

class MolDyn {

	public :

	MolDyn();
	MolDyn(bool, std::string old = std::string());
   	~MolDyn() {;} 

	int Get_steps() {return m_steps;}
	int Get_throws() {return m_throws;}
	int Get_blocks() {return m_blocks;}
	int Get_iprint() {return m_iprint;}

	void Move();
//	void ConfXYZ(int);
	void Measure();
	void Measure(double stime[4]);
	void ConfFinal(std::string, std::string);
	double Force(int, int);
	double Pbc(double);

	private :

	int m_blocks, m_steps, m_parts, m_iprint, m_seed, m_throws;
	double (*m_x), (*m_y), (*m_z), // x
		   (*m_xold), (*m_yold), (*m_zold), // xold
		   (*m_vx), (*m_vy), (*m_vz); // v
	double m_dt, m_rcut, m_box, m_rho, m_temp;

	std::ifstream openfile(std::string);

//	double acc,att;
//	int iv,ik,it,ie;
//	double energy;
//	double vol; penso non usato altrove
//	const int m_props;
};

#endif
