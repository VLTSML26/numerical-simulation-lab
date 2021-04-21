/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "MolDyn.h"

using namespace std;

MolDyn::MolDyn(bool equilibrate, string old) {

//	m_props = 4;
	m_seed = 1;
	m_blocks = 100;
	srand(m_seed);

	ifstream read = openfile("input.dat");
//	read.open("input.dat");
	read >> m_temp >> m_parts >> m_rho >> m_rcut >> m_dt >> m_steps >> m_iprint;
	read.close();

//	cout << m_temp << m_parts << m_rho << m_rcut << m_dt << m_steps << m_iprint << endl;;
	double vol = (double) m_parts / m_rho;
	m_box = pow(vol, 1./3.);
	m_throws = 100;

	/* per me questi sono inutili
	iv = 0;
  	ik = 1;
  	ie = 2;
  	it = 3;
	*/

	m_x = new double[m_parts];
	m_y = new double[m_parts];
	m_z = new double[m_parts];

	read = openfile("config.0");
  	for(int i=0; i<m_parts; ++i) {
    	read >> m_x[i] >> m_y[i] >> m_z[i];
    	m_x[i] *= m_box;
    	m_y[i] *= m_box;
    	m_z[i] *= m_box;
  	}
  	read.close();
	//cout << "MolDyn() finito" << endl;

	m_vx = new double[m_parts];
	m_vy = new double[m_parts];
	m_vz = new double[m_parts];

	m_xold = new double[m_parts];
	m_yold = new double[m_parts];
	m_zold = new double[m_parts];

	if(old.empty()) {
   		double sumv[3] = {0.};
    	for(int i=0; i<m_parts; ++i) {
			//cout << "old è empty" << endl;
     		m_vx[i] = rand()/double(RAND_MAX) - 0.5;
     		m_vy[i] = rand()/double(RAND_MAX) - 0.5;
     		m_vz[i] = rand()/double(RAND_MAX) - 0.5;

			//cout << "alora " << endl; //m_vx[i] << endl;
     		sumv[0] += m_vx[i];
     		sumv[1] += m_vy[i];
     		sumv[2] += m_vz[i];
   		}
	
		//cout << "ok" << endl;
		// per evitare drift del CM
   		for(int idim=0; idim<3; ++idim)
			sumv[idim] /= (double) m_parts;
   		double sumv2 = 0.;
   		for(int i=0; i<m_parts; ++i) {
     		m_vx[i] -= sumv[0];
     		m_vy[i] -= sumv[1];
     		m_vz[i] -= sumv[2];
	
     		sumv2 += m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i];
   		}
   		sumv2 /= (double) m_parts;
		
   		double fs = sqrt(3 * m_temp / sumv2); // fs = fattore di scala velocità
   		for(int i=0; i<m_parts; ++i) {
     		m_vx[i] *= fs;
     		m_vy[i] *= fs;
     		m_vz[i] *= fs;

			m_xold[i] = Pbc(m_x[i] - m_vx[i] * m_dt);
     		m_yold[i] = Pbc(m_y[i] - m_vy[i] * m_dt);
     		m_zold[i] = Pbc(m_z[i] - m_vz[i] * m_dt);
		}  	
	}

	else {
  		ifstream read = openfile("old.0");
		for(int j=0; j<m_parts; ++j) {
			read >> m_xold[j] >> m_yold[j] >> m_zold[j];
			m_xold[j] *= m_box; 
			m_yold[j] *= m_box;
			m_zold[j] *= m_box;
		}
		read.close();

		if(equilibrate) {
  			double xnew[m_parts], ynew[m_parts], znew[m_parts], 
				   fx[m_parts], fy[m_parts], fz[m_parts];

  			for(int i=0; i<m_parts; ++i) {
    			fx[i] = Force(i, 0);
    			fy[i] = Force(i, 1);
    			fz[i] = Force(i, 2);
  			}
		
  			for(int i=0; i<m_parts; ++i) { 
    			xnew[i] = Pbc(2. * m_x[i] - m_xold[i] + fx[i] * m_dt * m_dt);
    			ynew[i] = Pbc(2. * m_y[i] - m_yold[i] + fy[i] * m_dt * m_dt);
    			znew[i] = Pbc(2. * m_z[i] - m_zold[i] + fz[i] * m_dt * m_dt);
	
    			m_vx[i] = Pbc(xnew[i] - m_x[i]) / m_dt;
    			m_vy[i] = Pbc(ynew[i] - m_y[i]) / m_dt;
    			m_vz[i] = Pbc(znew[i] - m_z[i]) / m_dt;
			}
	
			double t = 0.;
  			for(int i=0; i<m_parts; ++i) 
				t += 0.5 * (m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i]);
			
    		double T = (2./3.) * t / (double) m_parts; // temperatura
			double fscale = T / m_temp;
			cout << fscale << endl;
	
			for(int i=0; i<m_parts; ++i) {
				m_vx[i] /= sqrt(fscale);
				m_vy[i] /= sqrt(fscale);
				m_vz[i] /= sqrt(fscale);
			}
	
  			for(int i=0; i<m_parts; ++i) {
    			m_xold[i] = Pbc(xnew[i] - m_dt * m_vx[i]);
    			m_yold[i] = Pbc(ynew[i] - m_dt * m_vy[i]);
    			m_zold[i] = Pbc(znew[i] - m_dt * m_vz[i]);
	
    			m_x[i] = xnew[i];
    			m_y[i] = ynew[i];
    			m_z[i] = znew[i];
			}
		}
	}

	//cout << "MolDyn(bool, string) finito" << endl;
}	

void MolDyn::Move() {
  	double xnew, ynew, znew, fx[m_parts], fy[m_parts], fz[m_parts];

  	for(int i=0; i<m_parts; ++i) {
    	fx[i] = Force(i, 0);
    	fy[i] = Force(i, 1);
    	fz[i] = Force(i, 2);
  	}

  	for(int i=0; i<m_parts; ++i) {
    	xnew = Pbc(2. * m_x[i] - m_xold[i] + fx[i] * m_dt * m_dt);
    	ynew = Pbc(2. * m_y[i] - m_yold[i] + fy[i] * m_dt * m_dt);
    	znew = Pbc(2. * m_z[i] - m_zold[i] + fz[i] * m_dt * m_dt);

    	m_vx[i] = Pbc(xnew - m_xold[i]) / (2. * m_dt);
    	m_vy[i] = Pbc(ynew - m_yold[i]) / (2. * m_dt);
    	m_vz[i] = Pbc(znew - m_zold[i]) / (2. * m_dt);

    	m_xold[i] = m_x[i];
    	m_yold[i] = m_y[i];
    	m_zold[i] = m_z[i];

    	m_x[i] = xnew;
    	m_y[i] = ynew;
    	m_z[i] = znew;
	}
  	
	//cout << "Move() finito" << endl;
	return;
}

void MolDyn::Measure() {
//	int bin;
  	double v = 0., t = 0.;
  	ofstream measure;

  	measure.open("measure.out", ios::app);

	// energia potenziale
  	for(int i=0; i<m_parts-1; ++i) {
    	for(int j=i+1; j<m_parts; ++j) {
     		double dx = Pbc(m_xold[i] - m_xold[j]); // here I use old configurations [old = r(t)]
     		double dy = Pbc(m_yold[i] - m_yold[j]); // to be compatible with EKin which uses v(t)
     		double dz = Pbc(m_zold[i] - m_zold[j]); // => EPot should be computed with r(t)

     		double dr = dx*dx + dy*dy + dz*dz;
     		dr = sqrt(dr);

     		if(dr < m_rcut) {
       			double vij = 4./pow(dr, 12) - 4./pow(dr, 6);
				v += vij;
			}
		}          
	}

	// energia cinetica
  	for(int i=0; i<m_parts; ++i) 
		t += 0.5 * (m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i]);
   
    double stima_epot = v / (double) m_parts;
    double stima_ekin = t / (double) m_parts; 
    double stima_temp = (2./3.) * t / (double) m_parts; 
    double stima_etot = (t + v) / (double) m_parts;
    	
	measure << stima_temp << "\t" << stima_etot << "\t"
			<< stima_ekin << "\t" << stima_epot << "\n";

	measure.close();
    return;
}

void MolDyn::Measure(double stime[4]) {
  	double v = 0., t = 0.;

  	for(int i=0; i<m_parts-1; ++i) {
    	for(int j=i+1; j<m_parts; ++j) {
     		double dx = Pbc(m_xold[i] - m_xold[j]);
     		double dy = Pbc(m_yold[i] - m_yold[j]);
     		double dz = Pbc(m_zold[i] - m_zold[j]);

     		double dr = dx*dx + dy*dy + dz*dz;
     		dr = sqrt(dr);

     		if(dr < m_rcut) {
       			double vij = 4./pow(dr, 12) - 4./pow(dr, 6);
				v += vij;
			}
		}          
	}

  	for(int i=0; i<m_parts; ++i) 
		t += 0.5 * (m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i]);
   
    double stima_epot = v / (double) m_parts;
    double stima_ekin = t / (double) m_parts; 
    double stima_temp = (2./3.) * t / (double) m_parts; 
    double stima_etot = (t + v) / (double) m_parts;
    	
	stime[0] = stima_temp;
	stime[1] = stima_etot;
	stime[2] = stima_ekin;
	stime[3] = stima_epot;

    return;
}

void MolDyn::ConfFinal(string conf, string old) {
	ofstream write;

  	write.open(conf);
  	for(int i=0; i<m_parts; ++i) {
    	write << m_x[i]/m_box << "\t" 
		      << m_y[i]/m_box << "\t" 
		      << m_z[i]/m_box << "\n";
  	}
  	write.close();

  	write.open(old);
  	for(int i=0; i<m_parts; ++i) {
    	write << m_xold[i]/m_box << "\t" 
	          << m_yold[i]/m_box << "\t" 
			  << m_zold[i]/m_box << "\n";
  	}
  	write.close();

  	return;
}

double MolDyn::Force(int ip, int idir) {
  	double f = 0.;
  	double dvec[3];

  	for(int i=0; i<m_parts; ++i) {
    	if(i != ip) {
      		dvec[0] = Pbc(m_x[ip] - m_x[i]);
      		dvec[1] = Pbc(m_y[ip] - m_y[i]);
      		dvec[2] = Pbc(m_z[ip] - m_z[i]);

      		double dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      		dr = sqrt(dr);

      		if(dr < m_rcut) {
        		f += dvec[idir] * (48./pow(dr, 14) - 24./pow(dr, 8)); // -Grad_ip V(r)
      		}
    	}
  	}
  
  	return f;
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

/*
MolDyn::MolDyn() {

//	m_props = 4;
	m_seed = 1;

	ifstream read = openfile("input.dat");
//	read.open("input.dat");
	read >> m_temp >> m_parts >> m_rho >> m_rcut >> m_dt >> m_steps >> m_iprint;
	read.close();

//	cout << m_temp << m_parts << m_rho << m_rcut << m_dt << m_steps << m_iprint << endl;;
	double vol = (double) m_parts / m_rho;
	m_box = pow(vol, 1./3.);

	*/
	/* per me questi sono inutili
	iv = 0;
  	ik = 1;
  	ie = 2;
  	it = 3;
	*/
/*
	m_x = new double[m_parts];
	m_y = new double[m_parts];
	m_z = new double[m_parts];

	read = openfile("config.0");
  	for(int i=0; i<m_parts; ++i) {
    	read >> m_x[i] >> m_y[i] >> m_z[i];
    	m_x[i] *= m_box;
    	m_y[i] *= m_box;
    	m_z[i] *= m_box;
  	}
  	read.close();
	cout << "MolDyn() finito" << endl;
}

MolDyn::MolDyn(bool equilibrate, string old) {
	MolDyn();
	srand(m_seed);

	m_vx = new double[m_parts];
	m_vy = new double[m_parts];
	m_vz = new double[m_parts];

	m_xold = new double[m_parts];
	m_yold = new double[m_parts];
	m_zold = new double[m_parts];

	if(old.empty()) {
   		double sumv[3] = {0.};
    	for(int i=0; i<m_parts; ++i) {
			cout << "old è empty" << endl;
     		m_vx[i] = rand()/double(RAND_MAX) - 0.5;
     		m_vy[i] = rand()/double(RAND_MAX) - 0.5;
     		m_vz[i] = rand()/double(RAND_MAX) - 0.5;

			cout << "alora " << endl; //m_vx[i] << endl;
     		sumv[0] += m_vx[i];
     		sumv[1] += m_vy[i];
     		sumv[2] += m_vz[i];
   		}
	
		cout << "ok" << endl;
		// per evitare drift del CM
   		for(int idim=0; idim<3; ++idim)
			sumv[idim] /= (double) m_parts;
   		double sumv2 = 0.;
   		for(int i=0; i<m_parts; ++i) {
     		m_vx[i] -= sumv[0];
     		m_vy[i] -= sumv[1];
     		m_vz[i] -= sumv[2];
	
     		sumv2 += m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i];
   		}
   		sumv2 /= (double) m_parts;
		
   		double fs = sqrt(3 * m_temp / sumv2); // fs = fattore di scala velocità
   		for(int i=0; i<m_parts; ++i) {
     		m_vx[i] *= fs;
     		m_vy[i] *= fs;
     		m_vz[i] *= fs;

			m_xold[i] = Pbc(m_x[i] - m_vx[i] * m_dt);
     		m_yold[i] = Pbc(m_y[i] - m_vy[i] * m_dt);
     		m_zold[i] = Pbc(m_z[i] - m_vz[i] * m_dt);
		}  	
	}

	else {
  		ifstream read = openfile("old.0");
		for(int j=0; j<m_parts; ++j) {
			read >> m_xold[j] >> m_yold[j] >> m_zold[j];
			m_xold[j] *= m_box; 
			m_yold[j] *= m_box;
			m_zold[j] *= m_box;
		}
		read.close();

		if(equilibrate) {
  			double xnew[m_parts], ynew[m_parts], znew[m_parts], 
				   fx[m_parts], fy[m_parts], fz[m_parts];

  			for(int i=0; i<m_parts; ++i) {
    			fx[i] = Force(i, 0);
    			fy[i] = Force(i, 1);
    			fz[i] = Force(i, 2);
  			}
		
  			for(int i=0; i<m_parts; ++i) { 
    			xnew[i] = Pbc(2. * m_x[i] - m_xold[i] + fx[i] * m_dt * m_dt);
    			ynew[i] = Pbc(2. * m_y[i] - m_yold[i] + fy[i] * m_dt * m_dt);
    			znew[i] = Pbc(2. * m_z[i] - m_zold[i] + fz[i] * m_dt * m_dt);
	
    			m_vx[i] = Pbc(xnew[i] - m_x[i]) / m_dt;
    			m_vy[i] = Pbc(ynew[i] - m_y[i]) / m_dt;
    			m_vz[i] = Pbc(znew[i] - m_z[i]) / m_dt;
			}
	
			double t = 0.;
  			for(int i=0; i<m_parts; ++i) 
				t += 0.5 * (m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i]);
			
    		double T = (2./3.) * t / (double) m_parts; // temperatura
			double fscale = T / m_temp;
			cout << fscale << endl;
	
			for(int i=0; i<m_parts; ++i) {
				m_vx[i] /= sqrt(fscale);
				m_vy[i] /= sqrt(fscale);
				m_vz[i] /= sqrt(fscale);
			}
	
  			for(int i=0; i<m_parts; ++i) {
    			m_xold[i] = Pbc(xnew[i] - m_dt * m_vx[i]);
    			m_yold[i] = Pbc(ynew[i] - m_dt * m_vy[i]);
    			m_zold[i] = Pbc(znew[i] - m_dt * m_vz[i]);
	
    			m_x[i] = xnew[i];
    			m_y[i] = ynew[i];
    			m_z[i] = znew[i];
			}
		}
	}

	cout << "MolDyn(bool, string) finito" << endl;
}*/	
