/*
 * LeeWiswall.hpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 */
/*
 * LeeWiswall.cpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */
///. Header and cpp merged by Ali Raeini

#include <mpi.h>
#include <iostream>
#include "string.h"
#include <algorithm>
#include <thread>

#ifndef DISTR_PAR_NELDERMEAD_HPP_
#define DISTR_PAR_NELDERMEAD_HPP_
#define SIMPLEX(i,j) simplex[((indices[(i)])*nx) + (j) ]




inline void onestep(int* feval, int* itrUpdate, double* xCur, const double xM[], double* fCur,const double fPrev, const double best, double (*obj_function)(double *, int), const int nx, double rho, double xi, double gamout, double gamain)
{

	double* xR = new double[nx];
	double* xE = new double[nx];
	double* xC = new double[nx];
	const double fCurr=*fCur;

	// compute reflection and store function value in fAR
	for (int i = 0; i < nx; i++)     xR[i] = (1 + rho) * xM[i] - rho * xCur[i];
	double fAR = obj_function(xR, nx);  (*feval)++;

	if(best <= fAR && fAR <= fPrev) {      // accept reflection point
		memmove(xCur, xR, nx * sizeof(double));
		*fCur = fAR;
	} else if(fAR < best) {
		for (int i = 0; i < nx; i++)    xE[i] = (1 + xi) * xM[i] - xi * xCur[i];  // test for expansion
		double fAE = obj_function(xE, nx);  (*feval)++;
		if(fAE < fAR) {                      // accept expansion point
			 memmove(xCur, xE, nx * sizeof(double));
			 *fCur = fAE;
		} else {                           // eventual accept reflection point
			 memmove(xCur, xR, nx * sizeof(double));
			 *fCur = fAR;
		}
	} else if(fPrev <=fAR && fAR < fCurr) {

		for (int i = 0; i < nx; i++)    xC[i] = (1 + gamout) * xM[i] - gamout * xCur[i];  // do outside contraction
		double fAC = obj_function(xC, nx);  (*feval)++;
		if(fAC <= fAR) {  // accept outside contraction point
			 memmove(xCur, xC, nx * sizeof(double));
			 *fCur = fAC;
		} else {   // shrink
			 if(fAR < fCurr) {
				  memmove(xCur, xR, nx * sizeof(double)); *itrUpdate=2; // just move the memory, do not update
				 *fCur = fAR;
			 }
		}
	} else {

		for (int i = 0; i < nx; i++)    xC[i] = (1 - gamain) * xM[i] + gamain * xCur[i];   // do inside contraction
		double fAC = obj_function(xC, nx);  (*feval)++;
		if(fAC < fCurr) {
			 // accept inside contraction point
			memmove(xCur, xC, nx * sizeof(double));
			*fCur = fAC;
		} else if(fAR < fCurr) { // shrink
				  memmove(xCur, xR, nx * sizeof(double)); *itrUpdate=2; // just move the memory, do not update
				  *fCur = fAR;
		} else *itrUpdate=0;
	}
	delete[] xR;
	delete[] xE;
	delete[] xC;

}


class IndexSorter {
public:
    IndexSorter(double *arg) :
    objVals(arg) {};

    bool operator()(int i, int j) {
        return objVals[i] < objVals[j];
    }
private:
    double *objVals;
};



class LeeWiswall {
public:
    /**
     * Given initial guess, a step, the nx of the simplex,
     * and a pointer to an objective function.
     *
     * simplex: Array of doubles size nx*(nx+1)
     * nxs: Dimension of each of the nx+1 vectors.
     * obj_function: Pointer to the objective function, takes as argument
     *              a vector and its length, should return a double.
     */
	LeeWiswall(double *guess, double step, int nX,  double (*objf)(double *, int), int rnk, int nprc,  int pPerItr)
	 : nx(nX), rank(rnk), nproc(nprc), nPerIter(pPerItr), obj_function(objf)
	{
		 //init(guess, step, nX, obj_function, rank, nproc, nPerIter);
		 // Determine how many points are on the given processor, and their global indices. Based off this index update with the provided step nproc.

		 // points per processor (this "rounds" down)
		 pNxp1 = (nx + 1) / nproc;

		 // assign remainder to ranks 0, 1, 2, ...
		 if ((nx + 1) % nproc > rank) {
			  pNxp1++;
		 }

		 int globalFirstIndex = rank * ((nx + 1) / nproc)
			  + std::min((nx + 1) % nproc, rank);

		 indices = new int[pNxp1];
		 for (int i = 0; i < pNxp1; i++)    indices[i] = i;

		 this->simplex = new double[nx * pNxp1];
		 for (int i = 0; i < pNxp1; i++) {
			  for (int j = 0; j < nx; j++) {
					SIMPLEX(i,j) = guess[j];
					if (globalFirstIndex + i == j + 1)      SIMPLEX(i,j) += step;
			  }
		 }

		 objVals = new double[pNxp1];

		 //std::cout<<" nPerIter "<<nPerIter<<std::endl;

		 rho = 1.;// RHO > 0
		 xi = 2.*rho; // XI  > max(RHO, 1)
		 gamout = 0.7;// 0 < GAM < 1
		 gamain = 0.7;// 0 < GAM < 1
		 sig = 0.5; // 0 < SIG < 1
		 feval = 0;
	}

    // Deletes user passed simplex as well as all allocated memory.
	~LeeWiswall() {
		 delete[] indices;
		 delete[] simplex;
		 delete[] objVals;

	}


   /**
     * Find the point which minimizes the objective function, and return
     * an array of nx doubles. User is responsible to free that memory.
     *
     * Will return answer if less than 1e-6, or if max_iterations > 0, then after
     * max_iterations, whichever comes first.
     */
   double* solve(int max_iterations) {

	double* xM = new double[nx];

	//Compute objective function for initial values
	for (int i = 0; i < pNxp1; i++) {
		objVals[i] = obj_function(&SIMPLEX(i, 0), nx);
		feval++;
	}
	sort_simplex(); //Sort the simplex
	MPI_Allreduce(&(objVals[indices[0]]), &best, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	int iter = 0;

	while (best > 1e-6 && (max_iterations <= 0 || iter * nproc < max_iterations)) {

        int updated = 0;

        // compute centroid
        //if (iter % nPerIter == 0) {
			 for (int i = 0; i < nx; i++)       xM[i] = 0.;
			 double sumW=0.;
			 for (int i = 0; i < pNxp1 - nPerIter; i++)
			 {
				 double w=4+2*nPerIter+pNxp1 -0.2*i;
				 sumW+=w;
				  for (int j = 0; j < nx; j++)  xM[j] += w*SIMPLEX(i,j);         //Divide after. Possible overflow for large obj function values!
			 }
			 //for (int i = 0; i < nx; i++)      xM[i] /= (nx + 1 - nproc * nPerIter); //Divide from earlier, then compute
			 //for (int i = 0; i < pNxp1- nPerIter; i++)
				  //for (int j = 0; j < nx; j++)  xM[j] += SIMPLEX(i,j);         //Divide after. Possible overflow for large obj function values!
			 //for (int i = 0; i < nx; i++)      xM[i] /= (pNxp1-nPerIter); //Divide from earlier, then compute

			 double SumWreduce; // Reduce xM to Mreduce
			 double *Mreduce = new double[nx]; // Reduce xM to Mreduce
			 MPI_Allreduce(xM, Mreduce, nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			 MPI_Allreduce(&sumW, &SumWreduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			 //memmove(xM, Mreduce, nx * sizeof(double));
			 for (int i = 0; i < nx; i++)      xM[i] = Mreduce[i]/SumWreduce; //Divide from earlier, then compute
 			 delete[] Mreduce;
			 //for (int i = 0; i < nx; i++)      xM[i] = xM[i]/sumW; //Divide from earlier, then compute
		//}

		{

			#ifdef MULTITHREAD
			 std::vector<std::thread> threads; threads.reserve(nPerIter);
			 std::vector<int> itrUpdates(nPerIter,1);
			 std::vector<int> nfItrs(nPerIter,0);
			for (int itrd=0;itrd<nPerIter; ++itrd)
			{
				int iCur = pNxp1 - ((iter+itrd) % nPerIter) - 1;
				threads.emplace_back(std::thread( onestep, &(nfItrs[itrd]), &(itrUpdates[itrd]), &SIMPLEX(iCur,0),  xM, &objVals[indices[iCur]], objVals[indices[iCur-1]], best,  obj_function, nx, rho, xi, gamout, gamain));
			}
			for (int itrd=0;itrd<nPerIter; ++itrd)
			{
				threads[itrd].join();
				updated=updated|(itrUpdates[itrd]==1);
				feval+=nfItrs[itrd];
			}

			#else
			//for (int itrd=0;itrd<nPerIter; ++itrd)
			{
				 for (int itrd=0;itrd<nPerIter; ++itrd) {
					int iCur = pNxp1 - ((iter+itrd) % nPerIter) - 1;
					int nfItr(0), itrUpdate(1);
					double* xCur=&SIMPLEX(iCur,0);
					double fPrev=objVals[indices[iCur - 1]];
					onestep(&nfItr, &itrUpdate, xCur,  xM, &objVals[indices[iCur]], fPrev, best,  obj_function, nx, rho, xi, gamout, gamain) ;
					updated=updated|(itrUpdate==1);
					feval+=nfItr;
				}
			 }
		 #endif
		}
				iter+=nPerIter;

		//if ((iter % nPerIter) == 0) {
			int global_updated = 0;
			MPI_Allreduce(&updated, &global_updated, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (!global_updated) { //not one processor had an update, minimize
				 double *global_bestPoint = new double[nx];
				 global_best(global_bestPoint);
				 for (int i = 0; i < pNxp1; i++)
					 for (int j = 0; j < nx; j++) SIMPLEX(i,j) = sig * SIMPLEX(i,j) + (1.-sig)*global_bestPoint[j];
				 delete[] global_bestPoint;                //Re-eval all of the points
				 for (int i = 0; i < pNxp1; i++) {
					  objVals[indices[i]] = obj_function(&SIMPLEX(i, 0), nx); feval++;
				 }
			}
			sort_simplex(); //Sort the simplex
			//Find the new best
			best = objVals[indices[0]];
			//Update global min on all processors
			MPI_Allreduce(&(objVals[indices[0]]), &best, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		//}
        /*if (iter * nproc % 500 == 0 && rank == 0) {
         std::cout << iter << " " " " << best << std::endl;
        }*/

   }

	int total_feval(0);
	MPI_Reduce(&feval, &total_feval, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) std::cout << "Total Iterations: " << iter << ",   Function Evaluations: " << total_feval << std::endl;
	double *answer = new double[nx];
	global_best(answer);
	delete[] xM;

	return answer;
}




    // Set rho, otherwise assumed to be RHO
    void set_rho(double rho);
    // Set xi, otherwise assumed XI
    void set_xi(double xi);
    // Set gam, otherwise assumed GAM
    void set_gam(double gam);
    // Set sig, otherwise assumed to be SIG
    void set_sig(double sig);
    // Set minimimum improvement to do restart after some number of iterations
    void setRestartCriterion(int iterations, double improvement);

    double bestV() const {return best;};

private:
    void init(double *guess, double step, int nx,
              double (*obj_function)(double *vector, int nx), int rank, int nproc,
              int pointsPerIter);

     void global_best(double *globest) {
		 struct { double val; int rank; } myBest, globestVal;
		 myBest.val = objVals[indices[0]];
		 myBest.rank = rank;
		 MPI_Allreduce(&myBest, &globestVal, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
		 if (rank == globestVal.rank)  memmove(globest, &SIMPLEX(0, 0), nx * sizeof(double));
		 MPI_Bcast(globest, nx, MPI_DOUBLE, globestVal.rank, MPI_COMM_WORLD);
	};
    void print_simplex() {// Debugging purposes
		 for (int i = 0; i < pNxp1; i++)
		 {
			  for (int j = 0; j < nx; j++)
				 std::cout << SIMPLEX(i,j) << " ";
			 std::cout << std::endl;
		 }
		std::cout << std::endl;
	};
	void sort_simplex() {    std::sort(indices, indices + pNxp1, IndexSorter(objVals));  };
	double *simplex ;
	double *objVals;
	double rho, xi, gamain, gamout, sig, best;
	int *indices;
	int nx, pNxp1;
	int rank, nproc, nPerIter;
	int feval;
	double (*obj_function)(double *vector, int nx);

};
#endif /* DISTR_PAR_NELDERMEAD_HPP_ */
