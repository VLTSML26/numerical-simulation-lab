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

#include <iostream>
#include <string>
#include <functional>
#include <fstream>

class Random {

	public :

  	Random();
	Random(std::string);
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
	void Spherical3D(double, double[]);
//	void Metropolis(double[], int, double, std::function<double(double[])>, std::string);
	template <class T> void Metropolis(T xn[], T y[], double delta, std::function<double(T[])> p, std::string s) {
		int dim = sizeof *y;
		if(s == "Uniform") {
			for(int i=0; i<dim; i++)
				y[i] = Rannyu(xn[i] - delta, xn[i] + delta);
		} else if(s == "Gauss") {
			for(int i=0; i<dim; i++)
				y[i] = Gauss(xn[i], delta);
		} else {
			std::cout << "For now, Metropolis can function only with <Uniform> and <Gauss> sampling" << std::endl;
			return;
		}
		
		double alpha = p(y) / p(xn);
		if(alpha >= Rannyu()) {
			m_counter ++;
			for(int i=0; i<dim; i++)
				xn[i] = y[i];
		}
	}
	template <class T> void Metropolis(T xn[], int dim, double delta, std::function<double(T[])> p, std::string s) {
		T *y = new T[dim];
		if(s == "Uniform") {
			for(int i=0; i<dim; i++)
				y[i] = Rannyu(xn[i] - delta, xn[i] + delta);
		} else if(s == "Gauss") {
			for(int i=0; i<dim; i++)
				y[i] = Gauss(xn[i], delta);
		} else {
			std::cout << "For now, Metropolis can function only with <Uniform> and <Gauss> sampling" << std::endl;
			return;
		}
		
		double alpha = p(y) / p(xn);
		if(alpha >= Rannyu()) {
			m_counter ++;
			for(int i=0; i<dim; i++)
				xn[i] = y[i];
		}
		delete y;
	}
	//void Tune(double[], int, double&, std::function<double(double[])>, std::string);
	template <class T> void Tune(T xn[], int dim, double& delta, std::function<double(T[])> p, std::string s) {
		T *x = new T[dim];
		bool half = true;
		std::cout << "============ TUNING METROPOLIS ============" << std::endl;
		do {
			for(int idir=0; idir<dim; idir++)
				x[idir] = xn[idir];
			m_counter = 0;
			for(int i=0; i<m_equilibrate; i++)
				Metropolis<T>(x, dim, delta, p, s);
			double rate = (double) m_counter / m_equilibrate;
			std::cout << "rate: " << rate << "\tdelta: " << delta << std::endl;
			if(rate > 0.51)
				delta += 0.001; // fine tuning, better start with 0.1 when no idea of right delta
			else if(rate < 0.49)
				delta -= 0.001;
			else half = false;
		} while(half);
		delete x;
	}
  	bool Bool();

	private :

  	int m_m1, m_m2, m_m3, m_m4,
		m_l1, m_l2, m_l3, m_l4,
		m_n1, m_n2, m_n3, m_n4;
	int m_counter, m_equilibrate;

	std::ifstream openfile(std::string);
};

#endif 
