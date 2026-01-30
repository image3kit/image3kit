/*
 * LeeWiswall.cpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */
#include "NelderMeadPar.h"
#include <mpi.h>
#include <iostream>
#include "string.h"
#include <algorithm>
 #include <thread>

LeeWiswall::LeeWiswall(double *guess, double step, int nX,  double (*objf)(double *, int), int rnk, int nprc,  int pPerItr)
 : nx(nX), rank(rnk), nproc(nprc), nPerIter(pPerItr), obj_function(objf)  {
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


LeeWiswall::~LeeWiswall() {
    delete[] indices;
    delete[] simplex;
    delete[] objVals;

}


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

double* LeeWiswall::solve(int max_iterations) {

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
			 for (int i = 0; i < pNxp1 - nPerIter; i++)  {
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
			#define MULTITHREAD

			#ifdef MULTITHREAD
			 std::vector<std::thread> threads; threads.reserve(nPerIter);
			 std::vector<int> itrUpdates(nPerIter,1);
			 std::vector<int> nfItrs(nPerIter,0);
			for (int itrd=0;itrd<nPerIter; ++itrd)  {
				int iCur = pNxp1 - ((iter+itrd) % nPerIter) - 1;
				threads.emplace_back(std::thread( onestep, &(nfItrs[itrd]), &(itrUpdates[itrd]), &SIMPLEX(iCur,0),  xM, &objVals[indices[iCur]], objVals[indices[iCur-1]], best,  obj_function, nx, rho, xi, gamout, gamain));
			}
			for (int itrd=0;itrd<nPerIter; ++itrd)  {
				threads[itrd].join();
				updated=updated|(itrUpdates[itrd]==1);
				feval+=nfItrs[itrd];
			}

			#else
			//for (int itrd=0;itrd<nPerIter; ++itrd)  {
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

    if (rank == 0) {
        std::cout << "Total Iterations: " << iter << ",   Function Evaluations: " << total_feval << std::endl;
    }
    double *answer = new double[nx];
    global_best(answer);
	delete[] xM;

    return answer;
}


void LeeWiswall::global_best(double *global_best) {
    struct {
        double val;
        int rank;
    } myBest, global_bestVal;

    myBest.val = objVals[indices[0]];
    myBest.rank = rank;
    MPI_Allreduce(&myBest, &global_bestVal, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if (rank == global_bestVal.rank) {
        memmove(global_best, &SIMPLEX(0, 0), nx * sizeof(double));
    }
    MPI_Bcast(global_best, nx, MPI_DOUBLE, global_bestVal.rank, MPI_COMM_WORLD);

}




// Debugging purposes
void LeeWiswall::print_simplex() {
    for (int i = 0; i < pNxp1; i++)
        for (int j = 0; j < nx; j++)
            std::cout << SIMPLEX(i,j) << " ";
        std::cout << std::endl;
    std::cout << std::endl;
}

void LeeWiswall::sort_simplex() {
    std::sort(indices, indices + pNxp1, IndexSorter(objVals));
}
