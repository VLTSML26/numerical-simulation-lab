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
#include "NVE.h"

using namespace std;

NVE::NVE(string folder, bool equilibrate, string old) : MolDyn(folder) {

	m_props = 4 + m_nbins;
	m_folder = folder;

	// 0 - Potential Energy
	// 1 - Kinectic Energy
	// 2 - Temperature
	// 3 - Total Energy
	// 4 - 1st bin g(r)
	// ...
	// 4 + m_nbin - last bin g(r)
	m_walker = new double[m_props];
	for(int i=0; i<m_props; ++i)
			m_walker[i] = 0.;

	m_vx = new double[m_npart];
	m_vy = new double[m_npart];
	m_vz = new double[m_npart];

	m_xold = new double[m_npart];
	m_yold = new double[m_npart];
	m_zold = new double[m_npart];

	if(old.empty()) {
   		double sumv[3] = {0.};
    	for(int i=0; i<m_npart; ++i) {
     		m_vx[i] = m_rnd.Rannyu() - 0.5;
     		m_vy[i] = m_rnd.Rannyu() - 0.5;
     		m_vz[i] = m_rnd.Rannyu() - 0.5;

     		sumv[0] += m_vx[i];
     		sumv[1] += m_vy[i];
     		sumv[2] += m_vz[i];
   		}
	
   		for(int idim=0; idim<3; ++idim)
			sumv[idim] /= (double) m_npart;
   		double sumv2 = 0.;
   		for(int i=0; i<m_npart; ++i) {
     		m_vx[i] -= sumv[0];
     		m_vy[i] -= sumv[1];
     		m_vz[i] -= sumv[2];
	
     		sumv2 += m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i];
   		}
   		sumv2 /= (double) m_npart;
		
   		double fs = sqrt(3 * m_temp / sumv2); // fs = fattore di scala velocità
   		for(int i=0; i<m_npart; ++i) {
     		m_vx[i] *= fs;
     		m_vy[i] *= fs;
     		m_vz[i] *= fs;

			m_xold[i] = Pbc(m_x[i] - m_vx[i] * m_dt);
     		m_yold[i] = Pbc(m_y[i] - m_vy[i] * m_dt);
     		m_zold[i] = Pbc(m_z[i] - m_vz[i] * m_dt);
		}  	
	}

	else {
  		ifstream read = openfile(m_folder + "/old.0");
		for(int j=0; j<m_npart; ++j) {
			read >> m_xold[j] >> m_yold[j] >> m_zold[j];
			m_xold[j] *= m_box; 
			m_yold[j] *= m_box;
			m_zold[j] *= m_box;
		}
		read.close();

		if(equilibrate) {
  			double xnew[m_npart], ynew[m_npart], znew[m_npart], 
				   fx[m_npart], fy[m_npart], fz[m_npart];

			ofstream scale;

  			for(int i=0; i<m_npart; ++i) {
    			fx[i] = Force(i, 0);
    			fy[i] = Force(i, 1);
    			fz[i] = Force(i, 2);
  			}
		
  			for(int i=0; i<m_npart; ++i) { 
    			xnew[i] = Pbc(2. * m_x[i] - m_xold[i] + fx[i] * m_dt * m_dt);
    			ynew[i] = Pbc(2. * m_y[i] - m_yold[i] + fy[i] * m_dt * m_dt);
    			znew[i] = Pbc(2. * m_z[i] - m_zold[i] + fz[i] * m_dt * m_dt);
	
    			m_vx[i] = Pbc(xnew[i] - m_x[i]) / m_dt;
    			m_vy[i] = Pbc(ynew[i] - m_y[i]) / m_dt;
    			m_vz[i] = Pbc(znew[i] - m_z[i]) / m_dt;
			}
	
			double t = 0.;
  			for(int i=0; i<m_npart; ++i) 
				t += 0.5 * (m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i]);
			
    		double T = (2./3.) * t / (double) m_npart; // temperatura
			double fscale = T / m_temp;
			scale.open(m_folder + "/scale_factors.out", ios::app);
			scale << fscale << endl;
			scale.close();
	
			for(int i=0; i<m_npart; ++i) {
				m_vx[i] /= sqrt(fscale);
				m_vy[i] /= sqrt(fscale);
				m_vz[i] /= sqrt(fscale);
			}
	
  			for(int i=0; i<m_npart; ++i) {
    			m_xold[i] = Pbc(xnew[i] - m_dt * m_vx[i]);
    			m_yold[i] = Pbc(ynew[i] - m_dt * m_vy[i]);
    			m_zold[i] = Pbc(znew[i] - m_dt * m_vz[i]);
	
    			m_x[i] = xnew[i];
    			m_y[i] = ynew[i];
    			m_z[i] = znew[i];
			}
		}
	}
}	

NVE::~NVE() {
	delete m_walker;
	delete m_xold;
	delete m_yold;
	delete m_zold;
	delete m_vx;
	delete m_vy;
	delete m_vz;
}

void NVE::Move() {
  	double xnew, ynew, znew, fx[m_npart], fy[m_npart], fz[m_npart];

  	for(int i=0; i<m_npart; ++i) {
    	fx[i] = Force(i, 0);
    	fy[i] = Force(i, 1);
    	fz[i] = Force(i, 2);
  	}

  	for(int i=0; i<m_npart; ++i) {
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
}

void NVE::Measure() {
  	double v = 0., t = 0.;

	for(int k=0; k<m_props; ++k)
		m_walker[k] = 0.;

  	for(int i=0; i<m_npart-1; ++i) {
    	for(int j=i+1; j<m_npart; ++j) {
     		double dx = Pbc(m_xold[i] - m_xold[j]);
     		double dy = Pbc(m_yold[i] - m_yold[j]);
     		double dz = Pbc(m_zold[i] - m_zold[j]);
     		double dr = dx*dx + dy*dy + dz*dz;
     		dr = sqrt(dr);
			// g(r)
			if(m_nbins > 0) {
				if(dr < m_box / 2.) {
					int bin = (int) (2 * dr * m_nbins / m_box);
					m_walker[4 + bin] += 2;
				}
			}
			// potential energy
     		if(dr < m_rcut) {
       			double vij = 4./pow(dr, 12) - 4./pow(dr, 6);
				v += vij;
			}
		}          
	}

	// kinectic energy
  	for(int i=0; i<m_npart; ++i) 
		t += 0.5 * (m_vx[i] * m_vx[i] + m_vy[i] * m_vy[i] + m_vz[i] * m_vz[i]);
   
    m_walker[0] = v / (double) m_npart;
    m_walker[1] = t / (double) m_npart; 
    m_walker[2] = (2./3.) * t / (double) m_npart; 
    m_walker[3] = (t + v) / (double) m_npart;
}

/*
 * This method is usually called after Measure() in order to divide the 
 * bin's occupation value by double div, whence a good estimation for g(r)
 */
void NVE::Gofr(double sum[]) {
	for(int i=0; i<m_nbins; i++) {
		double bin_size = (m_box / 2.0) / (double)m_nbins;
		double r = i * bin_size;
		double deltaV = 4. * M_PI/3. * (pow((r+bin_size), 3.) - pow(r, 3.));
		double div = m_rho * m_npart * deltaV * m_throws;
		sum[4 + i] /= div;
	}
}

void NVE::ConfFinal(string conf, string old) {
	ofstream write;

  	write.open(conf);
  	for(int i=0; i<m_npart; ++i) {
    	write << m_x[i]/m_box << "\t" 
		      << m_y[i]/m_box << "\t" 
		      << m_z[i]/m_box << "\n";
  	}
  	write.close();

  	write.open(old);
  	for(int i=0; i<m_npart; ++i) {
    	write << m_xold[i]/m_box << "\t" 
	          << m_yold[i]/m_box << "\t" 
			  << m_zold[i]/m_box << "\n";
  	}
  	write.close();

  	return;
}

double NVE::Force(int ip, int idir) {
  	double f = 0.;
  	double dvec[3];

  	for(int i=0; i<m_npart; ++i) {
    	if(i != ip) {
      		dvec[0] = Pbc(m_x[ip] - m_x[i]);
      		dvec[1] = Pbc(m_y[ip] - m_y[i]);
      		dvec[2] = Pbc(m_z[ip] - m_z[i]);

      		double dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      		dr = sqrt(dr);

      		if(dr < m_rcut) {
        		f += dvec[idir] * (48./pow(dr, 14) - 24./pow(dr, 8));
      		}
    	}
  	}
  
  	return f;
}
