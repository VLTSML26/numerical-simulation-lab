/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "MolDyn_NVE.h"

using namespace std;

double error(double, double, int);

int main() {

  	Input();
  	int nconf = 1;
	const int N_blocks = 100;
       	const int N_throws = N_steps / N_blocks;
	
	// N.B. b sta per Block Method
	const string filename[4] = {"ave_epot.out", "ave_etot.out", "ave_ekin.out", "ave_temp.out"};
	ofstream b_out[4];

	for(int i=0; i<4; i++){
		b_out[i].open(filename[i]);
	}
	double sum_prog[4] = {0.}, sum2_prog[4] = {0.};

       	// ==== NOTA BENE ====	
	// il metodo a blocchi potrebbe intaccare le altre funzioni del codice, controllare meglio (prima però pranzo)
	
	for(int i=0; i<N_blocks; i++) {
     		double sum[4] = {0.};
		for(int j=0; j<N_throws; j++) {
			Move();
//     			if(i%iprint == 0) cout << "Number of time-steps: " << i << endl;
     			//if(j%10 == 0) { // configurazioni estremamente correlate: non serve calcolare proprietà ogni step
        			Measure();
				sum[0] += stima_pot;
				sum[1] += stima_etot;
				sum[2] += stima_kin;
				sum[3] += stima_temp;
	//		        ConfXYZ(nconf);
        			nconf ++;
     			//}
		}
		for(int k=0; k<4; k++) {
			sum_prog[k] += sum[k] / N_throws;
			sum2_prog[k] += sum[k] / N_throws * sum[k] / N_throws;
			double av = sum_prog[k] / (i + 1);
			double av2 = sum2_prog[k] / (i + 1);
			b_out[k] << av << "\t" << error(av, av2, i) << endl;
		}
  	}	
	for(int i=0; i<4; i++){
		b_out[i].open(filename[i]);
	}
	
	ConfFinal(); // scrive la conf. finale: serve per poter ricominciare da questa
  	return 0;
}

double error(double sum, double sum2, int n) {
  	if (n == 0) return 0;
  	else return sqrt((sum2 - sum*sum) / n);
}

void Input(void) {

	ifstream ReadInput, ReadConf;
  	double ep, ek, pr, et, vir; // ex: e = energia, x = tipo (e.g. p = potenziale)

 // 	cout << "Classic Lennard-Jones fluid" << endl;
 // 	cout << "Molecular dynamics simulation in NVE ensemble" << endl << endl;
 // 	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
 // 	cout << "The program uses Lennard-Jones units" << endl;

  	seed = 1;    // Set seed for random numbers
  	srand(seed); // Initialize random number generator
  
  	ReadInput.open("input.dat"); 

  	ReadInput >> temp;

  	ReadInput >> N_part;
 // 	cout << "Number of particles = " << N_part << endl;

  	ReadInput >> rho;
//  	cout << "Density of particles = " << rho << endl;
  	vol = (double)N_part/rho;
 // 	cout << "Volume of the simulation box = " << vol << endl;
  	box = pow(vol, 1./3.);
 // 	cout << "Edge of the simulation box = " << box << endl;

  	ReadInput >> rcut;
  	ReadInput >> dt;
  	ReadInput >> N_steps;
  	ReadInput >> iprint;

//  	cout << "The program integrates Newton equations with the Verlet method" << endl;
//  	cout << "Time step = " << dt << endl;
//  	cout << "Number of steps = " << N_steps << endl << endl;
  	ReadInput.close();

	//Prepare array for measurements
  	iv = 0; //Potential energy
  	ik = 1; //Kinetic energy
  	ie = 2; //Total energy
  	it = 3; //Temperature
  	n_props = 4; //Number of observables

	//Read initial configuration
//  	cout << "Read initial configuration from file config.0 " << endl << endl;
  	ReadConf.open("config.0");
  	for(int i=0; i<N_part; ++i) {
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] *= box; // x[i] = x[i] * box
    		y[i] *= box;
    		z[i] *= box;
  	}
  	ReadConf.close();

	//Prepare initial velocities
//   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.};
   	for(int i=0; i<N_part; ++i) {
     		vx[i] = rand()/double(RAND_MAX) - 0.5;
     		vy[i] = rand()/double(RAND_MAX) - 0.5;
     		vz[i] = rand()/double(RAND_MAX) - 0.5;

     		sumv[0] += vx[i];
     		sumv[1] += vy[i];
     		sumv[2] += vz[i];
   	}
   	for(int idim=0; idim<3; ++idim) sumv[idim] /= (double)N_part;
   	double sumv2 = 0.;
   	for(int i=0; i<N_part; ++i) {
     		vx[i] -= sumv[0];
     		vy[i] -= sumv[1];
     		vz[i] -= sumv[2];

     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)N_part;
	
   	double fs = sqrt(3 * temp / sumv2);   // fs = fattore di scala velocità
   	for(int i=0; i<N_part; ++i) {
     		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;

     		xold[i] = Pbc(x[i] - vx[i] * dt);
     		yold[i] = Pbc(y[i] - vy[i] * dt);
     		zold[i] = Pbc(z[i] - vz[i] * dt);
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
