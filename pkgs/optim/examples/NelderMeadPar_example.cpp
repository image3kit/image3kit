/*
 * NelderMead_Driver.cpp
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 * use points_per_iter of 1/3 of total points for expensive objFunc
 * run mpi parallel, no threading, for large  problems with cheap objFunc
 */

#include <cmath>
#include <iostream>
#include <mpi.h>
#include "NelderMeadPar.h"
#include "ObjFunction.h"
using namespace std;
int main(int argc, char **argv) {
	int rank = 0, nproc = 1;
	int nx = 300, pointsPerIter = 1, max_iterations = -1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc <= 2) {
		std::cerr << "Error: incorrect usage ./execname <prob_size><points per iter> <max iterations (optional)>\n";
		exit(1);
	}
	if (argc > 3) {
		max_iterations = atoi(argv[3]);
	} else {
		max_iterations = -1;
	}

	nx = atoi(argv[1]);
	pointsPerIter = atoi(argv[2]);
	double *guess = new double[nx];
	for (int i = 0; i < nx; i++)	guess[i] = 30.;

	//LeeWiswall *solver = new LeeWiswall(guess, 10., nx, objFunction3, rank, nproc, pointsPerIter);
	LeeWiswall *solver = new LeeWiswall(guess, 3., nx, objFunctionWebPy, rank, nproc, pointsPerIter);

	double t1, t2;
	t1 = MPI_Wtime();
	double *answer = solver->solve(max_iterations);
	t2 = MPI_Wtime();
	if (rank == 0) {
		std::cout << "Elapsed time during solve: " << t2 - t1 << std::endl;
		//for (int i = 0; i < nx - 1; i++)	std::cout << answer[i] << ", ";		std::cout << answer[nx - 1] << std::endl;

			std::cout << answer[0] << "... "<<answer[nx - 1]<< std::endl;;

		 std::cout << solver->bestV() << std::endl;
   }
	delete solver;
	delete[] answer;
	MPI_Finalize();
	return 0;
}
