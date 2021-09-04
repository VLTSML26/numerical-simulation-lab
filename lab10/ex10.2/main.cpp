/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#define N_cities		32
#define N_individuals	5000
#define N_steps			20
#define N_generations	20

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include "TSP.h"
#include "mpi.h"

using namespace std;

vector<City> Set(ofstream&, Random, string, int);
void TSP(string, ofstream&, ofstream&, ofstream&, Random, double, double, double, int);
double Accumulate(vector<Individual>);

int main(int argc, char *argv[]) {

	Random rnd;
	MPI_Init(&argc, &argv);

	// get number of threads
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size != 4) {
		MPI_Finalize();
		cerr << "Error: must have 4 processes" << endl;
		return 1;
	}

	// get rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Random rnd{SEED_DIR "/Primes", "seed/seed." + to_string(rank)};

	// initialize ofstream for both circle and square TSP: the program prints
	// 1 - the cities positions (pos)
	// 2 - the best path between them (best)
	// 3 - the cost of the best individual of the final generation
  	ofstream pos[2];
	ofstream best[2];
	ofstream mean[2];
	string type[2] = {"circle", "square"};

	// set the probabilities for each process
	double exp[2] = {2.5, 2.};
	double pcross[2] = {0.69, 0.7};
	double pmutate[2] = {0.33, 0.1};

	// iterate over type = circle and type = square
	for(int i=0; i<2; ++i) {
		pos[i].open("data/" + type[i] + "_pos_" + to_string(rank) + ".out");
		best[i].open("data/" + type[i] + "_best_" + to_string(rank) + ".out");
		mean[i].open("data/" + type[i] + "_mean_" + to_string(rank) + ".out");

		// calls TSP algorithm
		TSP(type[i], pos[i], best[i], mean[i], rnd, exp[i], pcross[i], pmutate[i], rank);
		pos[i].close();
		best[i].close();
	}

	MPI_Finalize();

	return 0;
}

/* 
 * This function returns a vector<City> of randomly places
 * cities on a type = "circle" or inside a type = "square".
 */
vector<City> Set(ofstream &pos, Random rnd, string type, int rank) {

	// it is different from previous codes, since processes 1,2 and 3
	// do not allocate memory while they write cities[i].
	// if you use push_back() you will end up with a segmentation fault
	// see line 97 (if rank == 0)
	// line 101 i do not use push_back()
	vector<City> cities(N_cities);

	if(type == "circle") {
		MPI_Barrier(MPI_COMM_WORLD);

		// process 0 generates the m_city placed on a circle
		if(rank == 0) {
			for(int i=0; i<N_cities; ++i) {
				double theta = rnd.Rannyu(0, 2*M_PI);
				City c(cos(theta), sin(theta), i);
				cities[i] = c;
				pos	<< c.Get_x() << "\t" << c.Get_y() << endl;
			}
		}

		// broadcast the cities to all threads
		MPI_Bcast(cities.data(), cities.size() * sizeof(decltype(cities)::value_type), MPI_BYTE, 0, MPI_COMM_WORLD);
	}
	else if(type == "square") {
		MPI_Barrier(MPI_COMM_WORLD);

		// process 0 generates the m_city placed on a square
		if(rank == 0) {
			for(int i=0; i<N_cities; ++i) {
				double x = rnd.Rannyu(0, 1);
				double y = rnd.Rannyu(0, 1);
				City c(x, y, i);
				cities[i] = c;
				pos	<< c.Get_x() << "\t" << c.Get_y() << endl;
			}
		}

		// broadcast the cities to all threads
		MPI_Bcast(cities.data(), cities.size() * sizeof(decltype(cities)::value_type), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		cout << "Error: for now, let's place the cities on a <circle> or a <square>" << endl;
		exit(0);
	}

	return cities;
}

/*
 * This function implements the MPI Traveling Salesman Algorithm
 */
void TSP(string type, ofstream &pos, ofstream &best, ofstream &out, Random rnd, double exp, double pcross, double pmutate, int rank) {

	// set cities (type = circle or type = square)
	vector<City> cities = Set(pos, rnd, type, rank);

	// initialize vector of N_individuals individuals
	vector<Individual> ind;
	for(int i=0; i<N_individuals; ++i) {							
		Individual io(cities, rnd);
		ind.push_back(io);
	}

	// generate population
	Population p(ind, rnd, N_cities, exp, pcross, pmutate);

	out << p.m_ind[0].Cost() << '\t' << Accumulate(p.m_ind) << endl;
	for(int i=0; i<N_steps; i++) {

		// evolve the population for N_generations
		for(int j=0; j<N_generations; j++) {
			p.Evolve();
			out << p.m_ind[0].Cost() << '\t' << Accumulate(p.m_ind) << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int cont;

		// thread 0 selects a random continent to exchange the best with
		// the other two threads will exchange their bests
		if(rank == 0)
			cont = (int) rnd.Rannyu(1, 4);
		MPI_Bcast(&cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
		vector<City> send{p.m_ind[0].m_city};
		vector<City> recv{send.size()};
		size_t size = send.size() * sizeof(decltype(send)::value_type);

		// the process with lower rank does send-recv, the one with higher rank does recv-send.
		if(rank == 0) {	// 0: send-recv with cont
			MPI_Sendrecv(send.data(), size, MPI_BYTE, cont, 0, recv.data(), size, MPI_BYTE, cont, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);

		} else if(rank == cont) {	// cont: recv-send with 0
			MPI_Recv(recv.data(), size, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
			MPI_Send(send.data(), size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		} else if(rank == 1) {	// 1, but not cont: send-recv with 2 or 3
			if(cont == 2)
				MPI_Sendrecv(send.data(), size, MPI_BYTE, 3, 0, recv.data(), size, MPI_BYTE, 3, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
			else
				MPI_Sendrecv(send.data(), size, MPI_BYTE, 2, 0, recv.data(), size, MPI_BYTE, 2, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
		} else if(rank == 2) {	// 2, but not cont: send-recv with 3, or recv-send with 1
			if(cont == 1) {
				MPI_Sendrecv(send.data(), size, MPI_BYTE, 3, 0, recv.data(), size, MPI_BYTE, 3, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
			} else {
				MPI_Recv(recv.data(), size, MPI_BYTE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
				MPI_Send(send.data(), size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
			}
		} else {	// 3, but not cont: recv-send with 1 or 2
			if(cont == 1) {
				MPI_Recv(recv.data(), size, MPI_BYTE, 2, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
				MPI_Send(send.data(), size, MPI_BYTE, 2, 0, MPI_COMM_WORLD);
			} else {
				MPI_Recv(recv.data(), size, MPI_BYTE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
				MPI_Send(send.data(), size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		p.m_ind[0].m_city = recv;
		p.Sort();
	}

	// write the best path to the output file
	for(int i=0; i<N_cities; i++)
		best << p.m_ind[0].m_city[i].Get_label() << endl;
}

double Accumulate(vector<Individual> ind) {
	return accumulate(ind.begin(), ind.begin() + N_individuals / 2, ind[0].Cost(), [](const double x, const Individual &io) {return x + io.Cost();}) / (N_individuals / 2);
}
