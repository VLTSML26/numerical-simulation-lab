/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "MolDyn.h"

using namespace std;

int main(int argc, char* argv[]) {

  	if(argc<2) {
	  	cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
	  	return -1;
	}

	if(string(argv[1]) == "start") init(false);
	else if(string(argv[1]) == "equilibrate") init(true, "old.0");
	else if(string(argv[1]) == "measure") init(false, "old.0");
	else {
		cerr << "Usage: " << argv[0] << " start/equilibrate/measure\n";
		return 1;
	}

	int nconf = 1;
  	for(int i=1; i <= N_steps; ++i) {
     		Move();
     		if(i%iprint == 0) cout << i/iprint * 10 << "%\r";
		cout.flush();
     		if(i%10 == 0) {
        		Measure();
//		        ConfXYZ(nconf); 
        		nconf += 1;
     		}
		if(i == N_steps - 1) ConfOld();
  	}
  	
	cout << endl;
	ConfFinal();
  	return 0;
}

void init(bool equilibrate, string old) {

	ifstream ReadInput, ReadConf, ReadOld;
  	double ep, ek, pr, et, vir; // ex: e = energia, x = tipo (e.g. p = potenziale)

  	seed = 1; 
  	srand(seed);
  
  	ReadInput.open("input.dat"); 
  	ReadInput >> temp >> N_part >> rho >> rcut >> dt >> N_steps >> iprint;
  	ReadInput.close();

	vol = (double)N_part/rho;
	box = pow(vol, 1./3.);

  	iv = 0; //Potential energy
  	ik = 1; //Kinetic energy
  	ie = 2; //Total energy
  	it = 3; //Temperature
  	n_props = 4; //Number of observables

  	ReadConf.open("config.0");
  	for(int i=0; i<N_part; ++i) {
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] *= box;
    		y[i] *= box;
    		z[i] *= box;
  	}
  	ReadConf.close();

     	if(old.empty()) {
   		double sumv[3] = {0.};
    		for(int i=0; i<N_part; ++i) {
     			vx[i] = rand()/double(RAND_MAX) - 0.5;
     			vy[i] = rand()/double(RAND_MAX) - 0.5;
     			vz[i] = rand()/double(RAND_MAX) - 0.5;
	
     			sumv[0] += vx[i];
     			sumv[1] += vy[i];
     			sumv[2] += vz[i];
   		}
	
		// per evitare drift del CM
   		for(int idim=0; idim<3; ++idim) sumv[idim] /= (double)N_part;
   		double sumv2 = 0.;
   		for(int i=0; i<N_part; ++i) {
     			vx[i] -= sumv[0];
     			vy[i] -= sumv[1];
     			vz[i] -= sumv[2];
	
     			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   		}
   		sumv2 /= (double)N_part;
		
   		double fs = sqrt(3 * temp / sumv2); // fs = fattore di scala velocità
   		for(int i=0; i<N_part; ++i) {
     			vx[i] *= fs;
     			vy[i] *= fs;
     			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * dt);
     			yold[i] = Pbc(y[i] - vy[i] * dt);
     			zold[i] = Pbc(z[i] - vz[i] * dt);
		}  	
	}
	else {
  		ReadOld.open("old.0");
		for(int j=0; j<N_part; ++j) {
			ReadOld >> xold[j] >> yold[j] >> zold[j];
			xold[j] *= box; // x[i] = x[i] * box
			yold[j] *= box;
			zold[j] *= box;
		}
		ReadOld.close();

		if(equilibrate) {
  			double xnew[N_part], ynew[N_part], znew[N_part], fx[m_part], fy[m_part], fz[m_part];

  			for(int i=0; i<N_part; ++i) { //Force acting on particle i
    				fx[i] = Force(i, 0);
    				fy[i] = Force(i, 1);
    				fz[i] = Force(i, 2);
  			}
		
  			for(int i=0; i<N_part; ++i) { //Verlet integration scheme
		
    				xnew[i] = Pbc( 2. * x[i] - xold[i] + fx[i] * dt * dt );
    				ynew[i] = Pbc( 2. * y[i] - yold[i] + fy[i] * dt * dt );
    				znew[i] = Pbc( 2. * z[i] - zold[i] + fz[i] * dt * dt );
		
    				vx[i] = Pbc(xnew[i] - x[i]) / dt;
    				vy[i] = Pbc(ynew[i] - y[i]) / dt;
    				vz[i] = Pbc(znew[i] - z[i]) / dt;
		
			}
	
			double t = 0.;
			
  			for(int i=0; i<N_part; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
			
    			double T = (2./3.) * t/(double)N_part; //Temperature
			double fscale = T / temp;
			
			cout << fscale << endl;
	
			for(int i=0; i<N_part; ++i) {
				vx[i] /= sqrt(fscale);
				vy[i] /= sqrt(fscale);
				vz[i] /= sqrt(fscale);
			}
	
  			for(int i=0; i<N_part; ++i) {
    				xold[i] = Pbc(xnew[i] - dt*vx[i]);
    				yold[i] = Pbc(ynew[i] - dt*vy[i]);
    				zold[i] = Pbc(znew[i] - dt*vz[i]);
		
    				x[i] = xnew[i];
    				y[i] = ynew[i];
    				z[i] = znew[i];
			}  	
		}
	}
	
   	return;
}

void Move(void) { //Move particles with Verlet algorithm
  	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  	for(int i=0; i<N_part; ++i) { //Force acting on particle i
    		fx[i] = Force(i, 0);
    		fy[i] = Force(i, 1);
    		fz[i] = Force(i, 2);
  	}

  	for(int i=0; i<N_part; ++i) { //Verlet integration scheme
	
    		xnew = Pbc( 2. * x[i] - xold[i] + fx[i] * dt * dt );
    		ynew = Pbc( 2. * y[i] - yold[i] + fy[i] * dt * dt );
    		znew = Pbc( 2. * z[i] - zold[i] + fz[i] * dt * dt );

    		vx[i] = Pbc(xnew - xold[i])/(2. * dt);
    		vy[i] = Pbc(ynew - yold[i])/(2. * dt);
    		vz[i] = Pbc(znew - zold[i])/(2. * dt);

    		xold[i] = x[i];
    		yold[i] = y[i];
    		zold[i] = z[i];

    		x[i] = xnew;
    		y[i] = ynew;
    		z[i] = znew;
  	}
  	return;
}

double Force(int ip, int idir) { //Compute forces as -Grad_ip V(r)
  	double f = 0.;
  	double dvec[3], dr;

  	for(int i=0; i<N_part; ++i) {
    		if(i != ip){
      			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      			dvec[1] = Pbc( y[ip] - y[i] );
      			dvec[2] = Pbc( z[ip] - z[i] );

      			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      			dr = sqrt(dr);

      			if(dr < rcut) {
        			f += dvec[idir] * (48./pow(dr, 14) - 24./pow(dr, 8)); // -Grad_ip V(r)
      			}
    		}
  	}
  
  	return f;
}

void Measure() { //Properties measurement
  	int bin;
  	double v, t, vij;
  	double dx, dy, dz, dr;
  	ofstream Epot, Ekin, Etot, Temp;

  	Epot.open("output_epot.dat",ios::app);
  	Ekin.open("output_ekin.dat",ios::app);
  	Temp.open("output_temp.dat",ios::app);
  	Etot.open("output_etot.dat",ios::app);

  	v = 0.; //reset observables
  	t = 0.;

	//cycle over pairs of particles
  	for(int i=0; i<N_part-1; ++i) {
    		for(int j=i+1; j<N_part; ++j) {

     			dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     			dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     			dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     			dr = dx*dx + dy*dy + dz*dz;
     			dr = sqrt(dr);

     			if(dr < rcut) {
       				vij = 4./pow(dr, 12) - 4./pow(dr, 6);

				//Potential energy
       				v += vij;
     			}
    		}          
  	}

	//Kinetic energy
  	for(int i=0; i<N_part; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    	stima_pot = v/(double)N_part; //Potential energy per particle
    	stima_kin = t/(double)N_part; //Kinetic energy per particle
    	stima_temp = (2./3.) * t/(double)N_part; //Temperature
    	stima_etot = (t+v)/(double)N_part; //Total energy per particle
    	
	Epot << stima_pot  << endl;
    	Ekin << stima_kin  << endl;
    	Temp << stima_temp << endl;
    	Etot << stima_etot << endl;

    	Epot.close();
    	Ekin.close();
    	Temp.close();
    	Etot.close();

    	return;
}


void ConfFinal(void) { //Write final configuration
  	ofstream WriteConf;

//  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");

  	for(int i=0; i<N_part; ++i) {
    		WriteConf << x[i]/box << "\t" <<  y[i]/box << "\t" << z[i]/box << endl;
  	}
  	WriteConf.close();
  	return;
}

void ConfXYZ(int nconf) { //Write configuration in .xyz format
  	ofstream WriteXYZ;

  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	WriteXYZ << N_part << endl;
  	WriteXYZ << "This is only a comment!" << endl;
  	for(int i=0; i<N_part; ++i) {
    		WriteXYZ << "LJ\t" << Pbc(x[i]) << "\t" <<  Pbc(y[i]) << "\t" << Pbc(z[i]) << endl;
  	}
  	WriteXYZ.close();
}

double Pbc(double r) {  // periodic boundary conditions
    	return r - box * rint(r/box);
 }

void ConfOld(void) { //Write final configuration
  	ofstream WriteOld;

  	WriteOld.open("old.final");

  	for(int i=0; i<N_part; ++i) {
    		WriteOld << x[i]/box << "\t" <<  y[i]/box << "\t" << z[i]/box << endl;
  	}
  	WriteOld.close();
  	return;
}
