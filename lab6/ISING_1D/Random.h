/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef _Random_h_
#define _Random_h_

#include <string>
#include <fstream>

class Random {

	public :

  	Random();
  	~Random() {;}

	int Get_counter() {return m_counter;}
	void Set_counter(int counter) {m_counter = counter;}

  	void SetRandom(int*, int, int);
  	void SaveSeed();
  	double Rannyu(void);
  	double Rannyu(double, double);
  	double Gauss(double, double);
  	double Exp(double);
  	double Lorentz(double, double);
  	double Line();
	void Metropolis(double[], double, double (*)(double[]), std::string);
	void Tune(double[], double&, double (*)(double[]), std::string);
  	bool Bool();

	private :

  	int m_m1, m_m2, m_m3, m_m4,
		m_l1, m_l2, m_l3, m_l4,
		m_n1, m_n2, m_n3, m_n4;
	int m_counter, m_equilibrate;

	std::ifstream openfile(std::string);
};

#endif 
