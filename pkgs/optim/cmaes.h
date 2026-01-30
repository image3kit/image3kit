/**
 * @file cmaes.h
 * @author Nikolaus Hansen, ported to C++ by Alexander Fabisch
 *
 * \mainpage
 * CMA-ES for non-linear function minimization.
 *
 * Copyright of C implementation by Nikolaus Hansen (e-mail:
 * hansen .AT. bionik.tu-berlin.de, hansen .AT. lri.fr), ported to C++ by
 * <a href="mailto:afabisch@googlemail.com"> Alexander Fabisch</a>.
 *
 * \section lgpl License
 *
 * Copyright 1996, 2003, 2007, 2011 Nikolaus Hansen, Alexander Fabisch
 *
 * This file is part of CMA-ESpp.
 *
 * CMA-ESpp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMA-ESpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CMA-ESpp. If not, see <http://www.gnu.org/licenses/>.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * \section purpose General Purpose
 *
 * The CMA-ES (Evolution Strategy with Covariance Matrix Adaptation) is a
 * robust search/optimization method. The goal is to minimize a given
 * objective function, \f$f: R^N \rightarrow R\f$. The CMA-ES should be
 * applied, if e.g. BFGS and/or conjugate gradient methods fail due to a
 * rugged search landscape (e.g. discontinuities, outliers, noise, local
 * optima, etc.). Learning the covariance matrix in the CMA-ES is similar
 * to learning the inverse Hessian matrix in a quasi-Newton method. On
 * smooth landscapes the CMA-ES is roughly ten times slower than BFGS,
 * assuming derivatives are not directly available. For up to \f$N=10\f$
 * parameters the simplex direct search method (Nelder & Mead) is
 * sometimes faster, but less robust than CMA-ES.  On considerably hard
 * problems the search (a single run) is expected to take between
 * \f$100\cdot N\f$ and \f$300\cdot N^2\f$ function evaluations. But you
 * might be lucky...
 *
 * \section application Application Remark
 *
 * The adaptation of the covariance matrix (e.g. by the CMA) is
 * equivalent to a general linear transformation of the problem
 * variables. Nevertheless, every problem specific knowledge about the
 * best problem transformation should be exploited before starting the
 * search procedure and an appropriate a priori transformation should be
 * applied to the problem. In particular a decision should be taken
 * whether variables, which are positive by nature, should be taken in
 * the log scale. A hard lower variable bound can also be realized by
 * taking the square. All variables should be re-scaled such that they
 * "live" in a similar search range width (for example, but not
 * necessarily between zero and one), such that the initial standard
 * deviation can be chosen the same for all variables.
 *
 *
 * \section links Links
 *  - http://www.lri.fr/~hansen/cmaesintro.html
 *  - http://www.lri.fr/~hansen/publications.html
 *
 * \section tut Tutorial
 * 	- http://www.lri.fr/~hansen/cmatutorial.pdf
 *
 * \section references References
 *
 * - Hansen, N, and S. Kern (2004).  Evaluating the CMA Evolution
 *   Strategy on Multimodal Test Functions. In: Eighth International
 *   Conference on Parallel Problem Solving from Nature PPSN VIII,
 *   Proceedings, pp. 282-291, Berlin: Springer
 *
 * - Hansen, N., S.D. Müller and P. Koumoutsakos (2003): Reducing the
 *   Time Complexity of the Derandomized Evolution Strategy with
 *   Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
 *   11(1).
 *
 * - Hansen, N. and A. Ostermeier (2001). Completely Derandomized
 *   Self-Adaptation in Evolution Strategies. Evolutionary Computation,
 *   9(2), pp. 159-195.
 *
 * - Hansen, N. and A. Ostermeier (1996). Adapting arbitrary normal
 *   mutation distributions in evolution strategies: The covariance
 *   matrix adaptation. In Proceedings of the 1996 IEEE International
 *   Conference on Evolutionary Computation, pp. 312-317.
 */


//The code is based on Hansen's ANSI C implementation of CMA-ES:  http://www.lri.fr/~hansen/cmaes_inmatlab.html#C
//7. Make sure that the scale of all objective parameter components of the function is somewhat similar and sigma corresponds to about 1/4 of the respective search intervals.


#pragma once

//#include "parameters.h"
//#include "random.h"
//#include "timings.h"
//#include "utils.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>





#include <cmath>
#include <limits>
#include <ostream>
#include <iostream>
#include <stdexcept>
#include <string>




#include <ctime>
#include <cassert>

#include <algorithm>
#include <cmath>
#include <string>

template<typename T>
T square(T d){  return d*d;  }

template<typename T>
T maxElement(const T* rgd, int len){  return *std::max_element(rgd, rgd + len);  }

template<typename T>
T minElement(const T* rgd, int len){  return *std::min_element(rgd, rgd + len);  }

template<typename T>
int maxIndex(const T* rgd, int len){  return std::max_element(rgd, rgd + len) - rgd;  }

template<typename T>
int minIndex(const T* rgd, int len){  return std::min_element(rgd, rgd + len) - rgd;  }

/** sqrt(a^2 + b^2) numerically stable. */
template<typename T>
T myhypot(T a, T b)
{
  const T fabsa = std::fabs(a), fabsb = std::fabs(b);
  if(fabsa > fabsb)  { const T r = b / a;  return fabsa*std::sqrt(T(1)+r*r);  }
  else if(b != T(0)) { const T r = a / b;  return fabsb*std::sqrt(T(1)+r*r);  }
  else               {                              return T(0); }
}


/**
 * @class Random
 * A pseudo random number generator.
 */
template<typename T>
class Random
{
  // variables for uniform()
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int rgrand[32];
  // variables for gauss()
  bool stored;
  T hold;
public:
  /**
   * @param seed use clock if 0
   */
  Random(long unsigned seed = 0) : hold(0.)
  {
    stored = false;
    if(seed < 1)
    {
      long int t = 100*time(0) + clock();
      seed = (long unsigned) (t < 0 ? -t : t);
    }
    start(seed);
  }
  /**
   * @param seed 0 == 1
   */
  void start(long unsigned seed)
  {
    stored = false;
    startseed = seed;
    if(seed < 1) seed = 1;
    aktseed = seed;
    for(int i = 39; i >= 0; --i)
    {
      long tmp = aktseed / 127773;
      aktseed = 16807* (aktseed - tmp* 127773) - 2836* tmp;
      if(aktseed < 0) aktseed += 2147483647;
      if(i < 32) rgrand[i] = aktseed;
    }
    aktrand = rgrand[0];
  }
  /**
   * @return (0,1)-normally distributed random number
   */
  T gauss(void)
  {
    if(stored)
    {
      stored = false;
      return hold;
    }
    stored = true;
    T x1, x2, rquad;
    do {
      x1 = 2.*uniform() - 1.;
      x2 = 2.*uniform() - 1.;
      rquad = x1*x1 + x2*x2;
    } while(rquad >= 1 || rquad <= 0);
    const T fac = std::sqrt(T(-2)*std::log(rquad)/rquad);
    hold = fac*x1;
    return fac*x2;
  }
  /**
   * @return (0,1)-uniform distributed random number
   */
  T uniform(void)
  {
    long tmp = aktseed / 127773;
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if(aktseed < 0)
      aktseed += 2147483647;
    tmp = aktrand / 67108865;
    aktrand = rgrand[tmp];
    rgrand[tmp] = aktseed;
    return (T) aktrand / T(2.147483647e9);
  }
};


/**
 * @class cpuTiming
 * A class for time measurements of the eigen decomposition.
 * cpuTiming measures overall time and times between calls of tic and toc. For
 * small time spans (up to 1000 seconds) CPU time via clock() is used. For large
 * time spans the fall-back to elapsed time from time() is used.
 * timings_update() must be called often enough to prevent the fallback.
 */
class cpuTiming
{
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted;

  double lastdiff;
  double tictoczwischensumme;

public:
  double totaltime; //! zeroed by re-calling timings_start
  double totaltotaltime;
  double tictoctime;
  double lasttictoctime;

  cpuTiming()
  {
    totaltotaltime = 0;
    start();
  }

  void start()
  {
    totaltime = 0;
    tictoctime = 0;
    lasttictoctime = 0;
    istic = 0;
    lastclock = clock();
    lasttime = time(NULL);
    lastdiff = 0;
    tictoczwischensumme = 0;
    isstarted = 1;
  }

  /**
   * @return time between last call of timings_*() and now
   */
  double update()
  {
    double diffc, difft;
    clock_t lc = lastclock; // measure CPU in 1e-6s
    time_t lt = lasttime; // measure time in s

    assert(isstarted == 1 && "timings_started() must be called before using timings... functions");

    lastclock = clock(); // measures at most 2147 seconds, where 1s = 1e6 CLOCKS_PER_SEC
    lasttime = time(NULL);
    diffc = (double) (lastclock - lc) / CLOCKS_PER_SEC; // is presumably in [-21??, 21??]
    difft = difftime(lasttime, lt); // is presumably an integer
    lastdiff = difft; // on the "save" side
    // use diffc clock measurement if appropriate
    if(diffc > 0 && difft < 1000)
      lastdiff = diffc;
    assert(lastdiff >= 0 && "BUG in time measurement");
    totaltime += lastdiff;
    totaltotaltime += lastdiff;
    if(istic)
    {
      tictoczwischensumme += lastdiff;
      tictoctime += lastdiff;
    }
    return lastdiff;
  }

  void tic()
  {
    assert(!istic && "Timingic called twice without toc");
    update();
    istic = 1;
  }

  double toc()
  {
    assert(istic && "Timingoc called without tic");
    update();
    lasttictoctime = tictoczwischensumme;
    tictoczwischensumme = 0;
    istic = 0;
    return lasttictoctime;
  }
};




template<typename T> class CMAES;

/**
 * @class Parameters
 * Holds all parameters that can be adjusted by the user.
 */
template<typename T>
class Parameters
{
  friend class CMAES<T>;
public:

  /* Input parameters. */
  //! Problem dimension, must stay constant.
  int N;
  //! Initial search space vector.
  T* xstart;
  //! A typical value for a search space vector.
  T* typicalX;
  //! Indicates that the typical x is the initial point.
  bool typicalXcase;
  //! Initial standard deviations.
  T* rgInitialStds;
  T* rgDiffMinChange;

  /* Termination parameters. */
  //! Maximal number of objective function evaluations.
  T stopMaxFunEvals;
  T facmaxeval;
  //! Maximal number of iterations.
  T stopMaxIter;
  //! Minimal fitness value. Only activated if flg is true.
  struct { bool flg; T val; } stStopFitness;
  //! Minimal value difference.
  T stopTolFun;
  //! Minimal history value difference.
  T stopTolFunHist;
  //! Minimal search space step size.
  T stopTolX;
  //! Defines the maximal condition number.
  T stopTolUpXFactor;

  /* internal evolution strategy parameters */

  int lambda; //!< Population size. Number of samples per iteration, at least two, generally > 4.

  int mu;   //!< Number of individuals used to recompute the mean.
  T mucov;

  T mueff;   //!< Variance effective selection mass, should be lambda/4.

  T* weights;   //!< Weights used to recombinate the mean sum up to one.

  T damps;   //!< Damping parameter for step-size adaption, d = inifinity or 0 means adaption  is turned off, usually close to one.

  T cs;   //!< cs^-1 (approx. n/3) is the backward time horizon for the evolution path ps and larger than one.

  T ccumcov;

  T ccov;   //!< ccov^-1 (approx. n/4) is the backward time horizon for the evolution path  pc and larger than one.

  T diagonalCov;
  struct { T modulo; T maxtime; } updateCmode;
  T facupdateCmode;


  enum Weights  //!< Determines the method used to initialize the weights.
  {
    UNINITIALIZED_WEIGHTS, LINEAR_WEIGHTS, EQUAL_WEIGHTS, LOG_WEIGHTS
  } weightMode;

  std::string resumefile;  //!< File that contains an optimization state that should be resumed.


  bool logWarnings;  //!< Set to true to activate logging warnings.

  std::ostream& logStream;  //!< Output stream that is used to log warnings, usually std::cerr.


  Parameters()
      : N(-1),
        xstart(0),
        typicalX(0),
        typicalXcase(false),
        rgInitialStds(0),
        rgDiffMinChange(0),
        stopMaxFunEvals(-1),
        facmaxeval(1.),
        stopMaxIter(-1.),
        stopTolFun(1e-7),
        stopTolFunHist(1e-7),
        stopTolX(0), // 1e-11*insigma would also be reasonable
        stopTolUpXFactor(1e3),
        lambda(-1),
        mu(-1),
        mucov(-1),
        mueff(-1),
        weights(0),
        damps(-1),
        cs(-1),
        ccumcov(-1),
        ccov(-1),
        facupdateCmode(1),
        weightMode(UNINITIALIZED_WEIGHTS),
        resumefile(""),
        logWarnings(false),
        logStream(std::cerr)
  {
    stStopFitness.flg = false;
    stStopFitness.val = -std::numeric_limits<T>::max();
    updateCmode.modulo = -1;
    updateCmode.maxtime = -1;
  }

  Parameters(const Parameters& parameters)
  {
    assign(parameters);
  }

  ~Parameters()
  {
    if(xstart)
      delete[] xstart;
    if(typicalX)
      delete[] typicalX;
    if(rgInitialStds)
      delete[] rgInitialStds;
    if(rgDiffMinChange)
      delete[] rgDiffMinChange;
    if(weights)
      delete[] weights;
  }

  Parameters& operator=(const Parameters& parameters)
  {
    assign(parameters);
    return *this;
  }

  /**
   * @param dimension Dimension of the search space \f$N\f$. No default
   *                  available, must be defined here or you have to set the
   *                  member manually.
   * @param inxstart Initial point in search space \f$x_0\f$, default (NULL) is
   *                 \f$(0.5,\ldots,0.5)^T + N(0, initialStdDev^2) \in R^N\f$.
   *                 This must be an array of size \f$N\f$.
   * @param inrgsigma Coordinatewise initial standard deviation of the sample
   *                  distribution (\f$\sigma \cdot \sqrt{C_{ii}} =
   *                  initialStdDev[i]\f$). The expected initial distance
   *                  between initialX and the optimum per coordinate should be
   *                  roughly initialStdDev. The entries should not differ by
   *                  several orders of magnitude. Default (NULL) is
   *                  \f$(0.3,\ldots,0.3)^T \in R^N\f$. This must be an array of
   *                  size \f$N\f$.
   */
  void init(int dimension = 0, const T* inxstart = 0, const T* inrgsigma = 0)
  {
    if(logWarnings)
    {
      if(!(xstart || inxstart || typicalX))
        logStream << "Warning: initialX undefined. typicalX = 0.5....5." << std::endl;
      if(!(rgInitialStds || inrgsigma))
        logStream << "Warning: initialStandardDeviations undefined. 0.3....3." << std::endl;
    }

    if(dimension <= 0 && N <= 0)
      throw std::runtime_error("Problem dimension N undefined.");
    else if(dimension > 0)
      N = dimension;

    if(weightMode == UNINITIALIZED_WEIGHTS) weightMode = LOG_WEIGHTS;//;LINEAR_WEIGHTS

    diagonalCov = 0; // default is 0, but this might change in future

    if(!xstart)
    {
      xstart = new T[N];
      if(inxstart)
      {
        for(int i=0; i<N; ++i)       xstart[i] = inxstart[i];
      }
      else if(typicalX)
      {
        typicalXcase = true;
        for(int i=0; i<N; ++i)       xstart[i] = typicalX[i];
      }
      else
      {
        typicalXcase = true;
        for(int i=0; i<N; ++i)       xstart[i] = 0.5;
      }
    }

    if(!rgInitialStds)
    {
      rgInitialStds = new T[N];
      if(inrgsigma)
      {
        for(int i=0; i<N; ++i)      rgInitialStds[i] = inrgsigma[i];
      }
      else
      {
        for(int i=0; i<N; ++i)      rgInitialStds[i] = T(0.3);
      }
    }

    supplementDefaults();
  }

private:
  void assign(const Parameters& p)
  {
    N = p.N;

    if(xstart)
      delete[] xstart;
    if(p.xstart)
    {
      xstart = new T[N];
      for(int i=0; i<N; ++i)     xstart[i] = p.xstart[i];
    }

    if(typicalX)
      delete[] typicalX;
    if(p.typicalX)
    {
      typicalX = new T[N];
      for(int i=0; i<N; ++i)
        typicalX[i] = p.typicalX[i];
    }

    typicalXcase = p.typicalXcase;

    if(rgInitialStds)    delete[] rgInitialStds;
    if(p.rgInitialStds)
    {
      rgInitialStds = new T[N];
      for(int i=0; i<N; ++i)      rgInitialStds[i] = p.rgInitialStds[i];
    }

    if(rgDiffMinChange) delete[] rgDiffMinChange;
    if(p.rgDiffMinChange)
    {
      rgDiffMinChange = new T[N];
      for(int i=0; i<N; ++i)     rgDiffMinChange[i] = p.rgDiffMinChange[i];
    }

    stopMaxFunEvals = p.stopMaxFunEvals;
    facmaxeval = p.facmaxeval;
    stopMaxIter = p.stopMaxIter;

    stStopFitness.flg = p.stStopFitness.flg;
    stStopFitness.val = p.stStopFitness.val;

    stopTolFun = p.stopTolFun;
    stopTolFunHist = p.stopTolFunHist;
    stopTolX = p.stopTolX;
    stopTolUpXFactor = p.stopTolUpXFactor;

    lambda = p.lambda;
    mu = p.mu;
    mucov = p.mucov;
    mueff = p.mueff;

    if(weights)  delete[] weights;
    if(p.weights)
    {
      weights = new T[mu];
      for(int i=0; i<mu; ++i)        weights[i] = p.weights[i];
    }

    damps = p.damps;
    cs = p.cs;
    ccumcov = p.ccumcov;
    ccov = p.ccov;
    diagonalCov = p.diagonalCov;

    updateCmode.modulo = p.updateCmode.modulo;
    updateCmode.maxtime = p.updateCmode.maxtime;

    facupdateCmode = p.facupdateCmode;

    weightMode = p.weightMode;

    resumefile = p.resumefile;
  }

  //! Supplements default parameter values.
  void supplementDefaults()
  {
    if(lambda < 2)      lambda =  4 + (int) (3.*log(0.+ N));
    if(mu <= 0)         mu = lambda / 2;
    if(!weights)        setWeights(weightMode);

    if(cs > 0)               cs *= (mueff + 2.) / (N + mueff + 3.);
    if(cs <= 0 || cs >= 1)   cs  = (mueff + 2.) / (N + mueff + 3.);

    if(ccumcov <= 0 || ccumcov > 1)   ccumcov = 4. / (N + 4);   ///aq orig=4. ->3.

    if(mucov < 1)                  mucov = mueff;
    T t1 = 2. / ((N + 1.4142)*(N + 1.4142));
    T t2 = (2.* mueff - 1.) / ((N + 2.)*(N + 2.) + mueff);
    t2 = (t2 > 1) ? 1 : t2;
    t2 = (1. / mucov)* t1 + (1. - 1. / mucov)* t2;
    if(ccov >= 0)              ccov *= t2;
    if(ccov < 0 || ccov > 1)   ccov = t2;

    //if(diagonalCov < 0)
      diagonalCov = 2 + 5. * N / sqrt((double) lambda);  ///aq orig=100. ->50.

    if(stopMaxFunEvals <= 0)      stopMaxFunEvals = facmaxeval * 900 * (N + 3)*(N + 3);
    else                          stopMaxFunEvals *= facmaxeval;

    if(stopMaxIter <= 0)          stopMaxIter = ceil((double) (stopMaxFunEvals / lambda));

    if(damps < T(0))   damps = T(1);
    damps = damps
        * (T(1) + T(2)*std::max(T(0), std::sqrt((mueff - T(1)) / (N + T(1))) - T(1)))
        * (T) std::max(T(0.3), T(1) - // modify for short runs
          (T) N / (T(1e-6) + std::min(stopMaxIter, stopMaxFunEvals / lambda)))
        + cs;

    if(updateCmode.modulo < 0)    updateCmode.modulo = 1. / ccov / (double) N / 10.;
    updateCmode.modulo *= facupdateCmode;
    if(updateCmode.maxtime < 0)   updateCmode.maxtime = 0.20; // maximal 20% of CPU-time
  }

  /**
   * Initializes the offspring weights.
   */
  void setWeights(Weights mode)
  {
    if(weights)
      delete[] weights;
    weights = new T[mu];
    switch(mode)
    {
    case LINEAR_WEIGHTS:
      for(int i=0; i<mu; ++i) weights[i] = mu - i;
      break;
    case EQUAL_WEIGHTS:
      for(int i=0; i<mu; ++i) weights[i] = 1;
      break;
    case LOG_WEIGHTS:
    default:
      for(int i=0; i<mu; ++i) weights[i] = log(mu + 1.) - log(i + 1.);
      break;
    }

    // normalize weights vector and set mueff
    T s1 = 0, s2 = 0;
    for(int i=0; i<mu; ++i)
    {
      s1 += weights[i];
      s2 += weights[i]*weights[i];
    }
    mueff = s1*s1/s2;
    for(int i=0; i<mu; ++i)
      weights[i] /= s1;

    if(mu < 1 || mu > lambda || (mu == lambda && weights[0] == weights[mu - 1]))
      throw std::runtime_error("setWeights(): invalid setting of mu or lambda");
  }
};






/**
 * @class CMAES
 * Evolution Strategies with Covariance Matrix Adaptation. The public interface
 * of the optimization algorithm.
 */
template<typename T>
class CMAES
{
public:

  /**
   * Keys for get().
   */
  enum GetScalar
  {
    NoScalar = 0,
    AxisRatio = 1,
    Eval = 2, Evaluations = 2,
    FctValue = 3, FuncValue = 3, FunValue = 3, Fitness = 3,
    FBestEver = 4,
    Generation = 5, Iteration = 5,
    MaxEval = 6, MaxFunEvals = 6, StopMaxFunEvals = 6,
    MaxGen = 7, MaxIter = 7, StopMaxIter = 7,
    MaxAxisLength = 8,
    MinAxisLength = 9,
    MaxStdDev = 10,
    MinStdDev = 11,
    Dim = 12, Dimension = 12,
    Lambda = 13, SampleSize = 13, PopSize = 13,
    Sigma = 14
  };

  /**
   * Keys for getPtr().
   */
  enum GetVector
  {
    NoVector = 0,
    DiagC = 1,
    DiagD = 2,
    StdDev = 3,
    XBestEver = 4,
    XBest = 5,
    XMean = 6
  };

  /**
   * Keys for writeToFile().
   */
  enum WriteKey
  {
    WCNone = 0,
    WKResume = 1,
    WKXMean = 2,
    WKC = 4,
    WKAll = 8,
    WKFewInfo = 16,
    WKFew = 32,
    WKEval = 64,
    WKFitness = 128,
    WKFBestEver = 256,
    WKCGeneration = 512,
    WKSigma = 1024,
    WKLambda = 2048,
    WKB = 4096,
    WKXBest = 8192,
    WKClock = 16384,
    WKDim = 32768
  };

private:

  std::string version;  //!< Implementation version.

  Random<T> rand;  //!< Random number generator.

  Parameters<T> params;  //!< CMA-ES parameters.

  T sigma;  //!< Step size.

  T* xmean;  //!< Mean x vector, "parent".

  T* xBestEver;  //! Best sample ever.

  T** population;  //! x-vectors, lambda offspring.

  int* index;  //!< Sorting index of sample population.

  T* funcValueHistory;  //!< History of function values.

  T chiN;

  T** C;  //!< Lower triangular matrix: i>=j for C[i][j].

  T** B;  //!< Matrix with normalize eigenvectors in columns.

  T* rgD;  //!< Axis lengths.

  T* pc;  //!< Anisotropic evolution path (for covariance).

  T* ps;  //!< Isotropic evolution path (for step length).

  T* xold;  //!< Last mean.

  T* output;  //!< Output vector.

  T* BDz;  //!< B*D*z.

  T* tempRandom;  //!< Temporary (random) vector used in different places.

  T* functionValues;  //!< Objective function values of the population.

  T* publicFitness;  //!< Public objective function value array returned by init().


  T gen;  //!< Generation number.

  enum {INITIALIZED, SAMPLED, UPDATED} state;  //!< Algorithm state.


  // repeatedly used for output
  T maxdiagC;
  T mindiagC;
  T maxEW;
  T minEW;

  bool eigensysIsUptodate;
  bool doCheckEigen; //!< control via signals.par
  T genOfEigensysUpdate;

  T dMaxSignifKond;
  T dLastMinEWgroesserNull;

  bool isResumeDone;

  time_t printtime;
  time_t writetime; //!< ideally should keep track for each output file
  time_t firstwritetime;
  time_t firstprinttime;

  std::string stopMessage; //!< A message that contains all matched stop criteria.

  std::string getTimeStr(void)
  {
    time_t tm = time(0);
    std::string timeStr(ctime(&tm));
    return timeStr.substr(0, 24);
  }

  /**
   * Calculating eigenvalues and vectors.
   * @param rgtmp (input) N+1-dimensional vector for temporal use.
   * @param diag (output) N eigenvalues.
   * @param Q (output) Columns are normalized eigenvectors.
   */
  void eigen(T* diag, T** Q, T* rgtmp)
  {
    assert(rgtmp && "eigen(): input parameter rgtmp must be non-NULL");

    if(C != Q) // copy C to Q
    {
      for(int i=0; i<params.N; ++i)
        for(int j = 0; j <= i; ++j)
          Q[i][j] = Q[j][i] = C[i][j];
    }

    householder(Q, diag, rgtmp);
    ql(diag, rgtmp, Q);
  }

  /**
   * Exhaustive test of the output of the eigendecomposition, needs O(n^3)
   * operations writes to error file.
   * @return number of detected inaccuracies
   */
  int checkEigen(T* diag, T** Q)
  {
    // compute Q diag Q^T and Q Q^T to check
    int res = 0;
    for(int i=0; i<params.N; ++i)
      for(int j = 0; j<params.N; ++j) {
        T cc = 0., dd = 0.;
        for(int k = 0; k < params.N; ++k)
        {
          cc += diag[k]*Q[i][k]*Q[j][k];
          dd += Q[i][k]*Q[j][k];
        }
        // check here, is the normalization the right one?
        const bool cond1 = fabs(cc - C[i > j ? i : j][i > j ? j : i]) / sqrt(C[i][i]* C[j][j]) > T(1e-10);
        const bool cond2 = fabs(cc - C[i > j ? i : j][i > j ? j : i]) > T(3e-14);
        if(cond1 && cond2)
        {
          std::stringstream s;
          s << i << " " << j << ": " << cc << " " << C[i > j ? i : j][i > j ? j : i]
              << ", " << cc - C[i > j ? i : j][i > j ? j : i];
          if(params.logWarnings)
            params.logStream << "eigen(): imprecise result detected " << s.str()
                << std::endl;
          ++res;
        }
        if(std::fabs(dd - (i == j)) > T(1e-10))
        {
          std::stringstream s;
          s << i << " " << j << " " << dd;
          if(params.logWarnings)
            params.logStream << "eigen(): imprecise result detected (Q not orthog.)"
                << s.str() << std::endl;
          ++res;
        }
      }
    return res;
  }

  /**
   * Symmetric tridiagonal QL algorithm, iterative.
   * Computes the eigensystem from a tridiagonal matrix in roughtly 3N^3 operations
   * code adapted from Java JAMA package, function tql2.
   * @param d input: Diagonale of tridiagonal matrix. output: eigenvalues.
   * @param e input: [1..n-1], off-diagonal, output from Householder
   * @param V input: matrix output of Householder. output: basis of
   *          eigenvectors, according to d
   */
  void ql(T* d, T* e, T** V)
  {
    const int n = params.N;
    T f(0);
    T tst1(0);
    const T eps(2.22e-16); // 2.^-52. = 2.22e-16

    // shift input e
    T* ep1 = e;
    for(T *ep2 = e+1, *const end = e+n; ep2 != end; ep1++, ep2++)
      *ep1 = *ep2;
    *ep1 = T(0); // never changed again

    for(int l = 0; l < n; l++)
    {
      // find small subdiagonal element
      T& el = e[l];
      T& dl = d[l];
      const T smallSDElement = std::fabs(dl) + std::fabs(el);
      if(tst1 < smallSDElement)
        tst1 = smallSDElement;
      const T epsTst1 = eps*tst1;
      int m = l;
      while(m < n)
      {
        if(std::fabs(e[m]) <= epsTst1) break;
        m++;
      }

      // if m == l, d[l] is an eigenvalue, otherwise, iterate.
      if(m > l)
      {
        do {
          T h, g = dl;
          T& dl1r = d[l+1];
          T p = (dl1r - g) / (T(2)*el);
          T r = myhypot(p, T(1));

          // compute implicit shift
          if(p < 0) r = -r;
          const T pr = p+r;
          dl = el/pr;
          h = g - dl;
          const T dl1 = el*pr;
          dl1r = dl1;
          for(int i = l+2; i < n; ++i) d[i] -= h;
          f += h;

          // implicit QL transformation.
          p = d[m];
          T c(1);
          T c2(1);
          T c3(1);
          const T el1 = e[l+1];
          T s(0);
          T s2(0);
          for(int i = m-1; i >= l; i--)
          {
            c3 = c2;
            c2 = c;
            s2 = s;
            const T& ei = e[i];
            g = c*ei;
            h = c*p;
            r = myhypot(p, ei);
            e[i+1] = s*r;
            s = ei/r;
            c = p/r;
            const T& di = d[i];
            p = c*di - s*g;
            d[i+1] = h + s*(c*g + s*di);

            // accumulate transformation.
            for(int k = 0; k < n; k++)
            {
              T& Vki1 = V[k][i+1];
              h = Vki1;
              T& Vki = V[k][i];
              Vki1 = s*Vki + c*h;
              Vki *= c; Vki -= s*h;
            }
          }
          p = -s*s2*c3*el1*el/dl1;
          el = s*p;
          dl = c*p;
        } while(std::fabs(el) > epsTst1);
      }
      dl += f;
      el = 0.;
    }
  }

  /**
   * Householder transformation of a symmetric matrix V into tridiagonal form.
   * Code slightly adapted from the Java JAMA package, function private tred2().
   * @param V input: symmetric nxn-matrix. output: orthogonal transformation
   *          matrix: tridiag matrix == V* V_in* V^t.
   * @param d output: diagonal
   * @param e output: [0..n-1], off diagonal (elements 1..n-1)
   */
  void householder(T** V, T* d, T* e)
  {
    const int n = params.N;

    for(int j = 0; j<n; j++)
    {
      d[j] = V[n - 1][j];
    }

    // Householder reduction to tridiagonal form

    for(int i = n - 1; i > 0; i--)
    {
      // scale to avoid under/overflow
      T scale = 0.;
      T h = 0.;
      for(T *pd = d, *const dend = d+i; pd != dend; pd++)
      {
        scale += std::fabs(*pd);
      }
      if(scale == 0.)
      {
        e[i] = d[i-1];
        for(int j = 0; j<i; j++)
        {
          d[j] = V[i-1][j];
          V[i][j] = 0.;
          V[j][i] = 0.;
        }
      }
      else
      {
        // generate Householder vector
        for(T *pd = d, *const dend = d+i; pd != dend; pd++)
        {
          *pd /= scale;
          h += *pd * *pd;
        }
        T& dim1 = d[i-1];
        T f = dim1;
        T g = f > 0 ? -std::sqrt(h) : std::sqrt(h);
        e[i] = scale*g;
        h = h - f* g;
        dim1 = f - g;
        memset((void *) e, 0, (size_t)i*sizeof(T));

        // apply similarity transformation to remaining columns
        for(int j = 0; j<i; j++)
        {
          f = d[j];
          V[j][i] = f;
          T& ej = e[j];
          g = ej + V[j][j]* f;
          for(int k = j + 1; k <= i - 1; k++)
          {
            T& Vkj = V[k][j];
            g += Vkj*d[k];
            e[k] += Vkj*f;
          }
          ej = g;
        }
        f = 0.;
        for(int j = 0; j<i; j++)
        {
          T& ej = e[j];
          ej /= h;
          f += ej* d[j];
        }
        T hh = f / (h + h);
        for(int j = 0; j<i; j++)
        {
          e[j] -= hh*d[j];
        }
        for(int j = 0; j<i; j++)
        {
          T& dj = d[j];
          f = dj;
          g = e[j];
          for(int k = j; k <= i - 1; k++)
          {
            V[k][j] -= f*e[k] + g*d[k];
          }
          dj = V[i-1][j];
          V[i][j] = 0.;
        }
      }
      d[i] = h;
    }

    // accumulate transformations
    const int nm1 = n-1;
    for(int i=0; i<nm1; ++i)
    {
      T h;
      T& Vii = V[i][i];
      V[n-1][i] = Vii;
      Vii = 1.;
      h = d[i+1];
      if(h != 0.)
      {
        for(int k = 0; k <= i; k++)
        {
          d[k] = V[k][i+1] / h;
        }
        for(int j = 0; j <= i; j++) {
          T g = 0.;
          for(int k = 0; k <= i; k++)
          {
            T* Vk = V[k];
            g += Vk[i+1]* Vk[j];
          }
          for(int k = 0; k <= i; k++)
          {
            V[k][j] -= g*d[k];
          }
        }
      }
      for(int k = 0; k <= i; k++)
      {
        V[k][i+1] = 0.;
      }
    }
    for(int j = 0; j<n; j++)
    {
      T& Vnm1j = V[n-1][j];
      d[j] = Vnm1j;
      Vnm1j = 0.;
    }
    V[n-1][n-1] = 1.;
    e[0] = 0.;
  }

  /**
   * Dirty index sort.
   */
  void sortIndex(const T* rgFunVal, int* iindex, int n)
  {
    int i, j;
    for(i = 1, iindex[0] = 0; i<n; ++i)
    {
      for(j = i; j > 0; --j)
      {
        if(rgFunVal[iindex[j - 1]] < rgFunVal[i])
          break;
        iindex[j] = iindex[j - 1]; // shift up
      }
      iindex[j] = i;
    }
  }

  void adaptC2(const int hsig)
  {
    const int N = params.N;
    bool diag = params.diagonalCov == 1 || params.diagonalCov >= gen;

    if(params.ccov != T(0))
    {
      // definitions for speeding up inner-most loop
      const T mucovinv = T(1)/params.mucov;
      const T commonFactor = params.ccov * (diag ? (N + T(1.5)) / T(3) : T(1));
      const T ccov1 = std::min(commonFactor*mucovinv, T(1));
      const T ccovmu = std::min(commonFactor*(T(1)-mucovinv), T(1)-ccov1);
      const T sigmasquare = sigma*sigma;
      const T onemccov1ccovmu = T(1)-ccov1-ccovmu;
      const T longFactor = (T(1)-hsig)*params.ccumcov*(T(2)-params.ccumcov);

      eigensysIsUptodate = false;

      // update covariance matrix
      for(int i=0; i<N; ++i)
        for(int j = diag ? i : 0; j <= i; ++j)
        {
          T& Cij = C[i][j];
          Cij = onemccov1ccovmu*Cij + ccov1 * (pc[i]*pc[j] + longFactor*Cij);
          for(int k = 0; k < params.mu; ++k)
          { // additional rank mu update
            const T* rgrgxindexk = population[index[k]];
            Cij += ccovmu*params.weights[k] * (rgrgxindexk[i] - xold[i])
                * (rgrgxindexk[j] - xold[j]) / sigmasquare;
          }
        }
      // update maximal and minimal diagonal value
      maxdiagC = mindiagC = C[0][0];
      for(int i=1; i < N; ++i)
      {
        const T& Cii = C[i][i];
        if(maxdiagC < Cii)
          maxdiagC = Cii;
        else if(mindiagC > Cii)
          mindiagC = Cii;
      }
    }
  }

  /**
   * Treats minimal standard deviations and numeric problems. Increases sigma.
   */
  void testMinStdDevs(void)
  {
    if(!this->params.rgDiffMinChange)
      return;

    for(int i=0; i<params.N; ++i)
      while(this->sigma*std::sqrt(this->C[i][i]) < this->params.rgDiffMinChange[i])
        this->sigma *= std::exp(T(0.05) + this->params.cs / this->params.damps);
  }

  /**
   * Adds the mutation sigma*B*(D*z).
   * @param x Search space vector.
   * @param eps Mutation factor.
   */
  void addMutation(T* x, T eps = 1.)
  {
    for(int i=0; i<params.N; ++i)
      tempRandom[i] = rgD[i]*rand.gauss();
    for(int i=0; i<params.N; ++i)
    {
      T sum = 0.;
      for(int j = 0; j<params.N; ++j)
        sum += B[i][j]*tempRandom[j];
      x[i] = xmean[i] + eps*sigma*sum;
    }
  }

  /**
   * This hack reads key words from input key for data to be written to
   * a file, see file signals.par as input file. The length of the keys
   * is mostly fixed. If the key phrase does not match the expectation the
   * output might be strange.
   */
  void writeToStream(int key, std::ostream& file)
  {
    if(key & WKResume)
    {
      file << std::endl << "# resume " << params.N << std::endl;
      file << "xmean" << std::endl;
      writeToStream(WKXMean, file);
      file << "path for sigma" << std::endl;
      for(int i=0; i<params.N; ++i)
        file << ps[i] << (i == params.N-1 ? "\n" : "\t");
      file << "path for C" << std::endl;
      for(int i=0; i<params.N; ++i)
        file << pc[i] << (i == params.N-1 ? "\n" : "\t");
      file << "sigma " << sigma << std::endl;
      // note than B and D might not be up-to-date
      file << "covariance matrix" << std::endl;
      writeToStream(WKC, file);
    }
    if(key & WKXMean)
    {
      for(int i=0; i<params.N; ++i)
        file << (i == 0 ? "" : "\t") << xmean[i];
      file << std::endl;
    }
    if(key & WKC)
    {
      for(int i=0; i<params.N; ++i)
        for(int j = 0; j <= i; ++j)
        {
          file << C[i][j];
          if(j == i)
            file << std::endl;
          else
            file << '\t';
        }
      file << std::endl;
    }
    if(key & WKAll)
    {
      time_t ti = time(0);
      file << std::endl << "# --------- " << asctime(localtime(&ti)) << std::endl;
      file << " N " << params.N << std::endl;
      file << "function evaluations " << (long) countevals << std::endl;
      file << "elapsed (CPU) time [s] " << std::setprecision(2) << eigenTimings.totaltotaltime << std::endl;
      file << "function value f(x)=" << population[index[0]][params.N] << std::endl;
      file << "maximal standard deviation " << sigma*std::sqrt(maxdiagC) << std::endl;
      file << "minimal standard deviation " << sigma*std::sqrt(mindiagC) << std::endl;
      file << "sigma " << sigma << std::endl;
      file << "axisratio " << (maxElement(rgD, params.N) / minElement(rgD, params.N)) << std::endl;
      file << "xbestever found after " << std::setprecision(0) << xBestEver[params.N+1]
          << "evaluations, function value " << xBestEver[params.N] << std::endl;
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << xBestEver[i] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "xbest (of last generation, function value " << population[index[0]][params.N] << ")" << std::endl;
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << population[index[0]][i] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "xmean" << std::endl;
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << xmean[i] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "Standard deviation of coordinate axes (sigma*sqrt(diag(C)))" << std::endl;
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << sigma*std::sqrt(C[i][i]) << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "Main axis lengths of mutation ellipsoid (sigma*diag(D))" << std::endl;
      for(int i=0; i<params.N; ++i)
          tempRandom[i] = rgD[i];
      std::sort(tempRandom, tempRandom + params.N);
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << sigma*tempRandom[params.N-1-i] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "Longest axis (b_i where d_ii=max(diag(D))" << std::endl;
      int k = maxIndex(rgD, params.N);
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << B[i][k] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << "Shortest axis (b_i where d_ii=max(diag(D))" << std::endl;
      k = minIndex(rgD, params.N);
      for(int i=0; i<params.N; ++i)
        file << " " << std::setw(12) << B[i][k] << (i % 5 == 4 || i == params.N-1 ? '\n' : ' ');
      file << std::endl;
    }
    if(key & WKFewInfo)
    {
      file << " Iter\tFevals\tFunction Value\tSigma\tMaxCoorDev\tMinCoorDev\t"
          << "AxisRatio\tMinDii\tTime in eig" << std::endl;
      file << std::endl;
    }
    if(key & WKFew)
    {
      file << (int) gen << "\t" << (int) countevals << "\t"
          << functionValues[index[0]] << "\t\t" << sigma << "  "
          << sigma*std::sqrt(maxdiagC) << "\t" << sigma*std::sqrt(mindiagC)
          << "\t" << std::scientific << std::setprecision(2)
          << std::sqrt(maxEW / minEW) << "\t" << std::sqrt(minEW)
          << "  " << eigenTimings.totaltotaltime;
      file << std::endl;
    }
    if(key & WKEval)
    {
      file << countevals;
      file << std::endl;
    }
    if(key & WKFitness)
    {
      for(int i=0; i<params.N; ++i)
        file << (i == 0 ? "" : "\t") << functionValues[index[i]];
      file << std::endl;
    }
    if(key & WKFBestEver)
    {
      file << xBestEver[params.N] << std::endl;
    }
    if(key & WKCGeneration)
    {
      file << gen << std::endl;
    }
    if(key & WKSigma)
    {
      file << sigma << std::endl;
    }
    if(key & WKLambda)
    {
      file << params.lambda << std::endl;
    }
    if(key & WKB)
    {
      int* iindex = new int[params.N];
      sortIndex(rgD, iindex, params.N);
      for(int i=0; i<params.N; ++i)
        for(int j = 0; j<params.N; ++j)
        {
          file << B[j][iindex[params.N-1-i]];
          if(j != params.N-1)
            file << '\t';
          else
            file << std::endl;
        }
      delete[] iindex;
      iindex = 0;
      file << std::endl;
    }
    if(key & WKXBest)
    {
      for(int i=0; i<params.N; ++i)
        file << (i == 0 ? "" : "\t") << population[index[0]][i];
      file << std::endl;
    }
    if(key & WKClock)
    {
      eigenTimings.update();
      file << eigenTimings.totaltotaltime << " " << eigenTimings.tictoctime
          << std::endl;
    }
    if(key & WKDim)
    {
      file << params.N;
      file << std::endl;
    }
  }

public:

  T countevals; //!< objective function evaluations
  cpuTiming eigenTimings;

  CMAES()
    : version("1.0alpha")
  {
  }

  /**
   * Releases the dynamically allocated memory, including that of the return
   * value of init().
   */
  ~CMAES()
  {
    delete[] pc;
    delete[] ps;
    delete[] tempRandom;
    delete[] BDz;
    delete[] --xmean;
    delete[] --xold;
    delete[] --xBestEver;
    delete[] --output;
    delete[] rgD;
    for(int i=0; i<params.N; ++i)
    {
      delete[] C[i];
      delete[] B[i];
    }
    for(int i=0; i<params.lambda; ++i)
      delete[] --population[i];
    delete[] population;
    delete[] C;
    delete[] B;
    delete[] index;
    delete[] publicFitness;
    delete[] --functionValues;
    delete[] --funcValueHistory;
  }

  /**
   * Initializes the CMA-ES algorithm.
   * @param parameters The CMA-ES parameters.
   * @return Array of size lambda that can be used to assign fitness values and
   *         pass them to updateDistribution(). Not that after the desctructor
   *         was called, the array is deleted.
   */
  T* init(const Parameters<T>& parameters)
  {
    params = parameters;

    stopMessage = "";

    T trace(0);
    for(int i=0; i<params.N; ++i)
      trace += params.rgInitialStds[i]*params.rgInitialStds[i];
    sigma = std::sqrt(trace/params.N);

    chiN = std::sqrt((T) params.N) * (T(1) - T(1)/(T(4)*params.N) + T(1)/(T(21)*params.N*params.N));
    eigensysIsUptodate = true;
    doCheckEigen = false;
    genOfEigensysUpdate = 0;
    isResumeDone = false;

    T dtest;
    for(dtest = T(1); dtest && dtest < T(1.1)*dtest; dtest *= T(2))
      if(dtest == dtest + T(1))
        break;
    dMaxSignifKond = dtest / T(1000); // not sure whether this is really save, 100 does not work well enough

    gen = 0;
    countevals = 0;
    state = INITIALIZED;
    dLastMinEWgroesserNull = T(1);
    printtime = writetime = firstwritetime = firstprinttime = 0;

    pc = new T[params.N];
    ps = new T[params.N];
    tempRandom = new T[params.N+1];
    BDz = new T[params.N];
    xmean = new T[params.N+2];
    xmean[0] = params.N;
    ++xmean;
    xold = new T[params.N+2];
    xold[0] = params.N;
    ++xold;
    xBestEver = new T[params.N+3];
    xBestEver[0] = params.N;
    ++xBestEver;
    xBestEver[params.N] = std::numeric_limits<T>::max();
    output = new T[params.N+2];
    output[0] = params.N;
    ++output;
    rgD = new T[params.N];
    C = new T*[params.N];
    B = new T*[params.N];
    publicFitness = new T[params.lambda];
    functionValues = new T[params.lambda+1];
    functionValues[0] = params.lambda;
    ++functionValues;
    const int historySize = 10 + (int) ceil(3.*10.*params.N/params.lambda);
    funcValueHistory = new T[historySize + 1];
    funcValueHistory[0] = (T) historySize;
    funcValueHistory++;

    for(int i=0; i<params.N; ++i)
    {
      C[i] = new T[i+1];
      B[i] = new T[params.N];
    }
    index = new int[params.lambda];
    for(int i=0; i<params.lambda; ++i)
        index[i] = i;
    population = new T*[params.lambda];
    for(int i=0; i<params.lambda; ++i)
    {
      population[i] = new T[params.N+2];
      population[i][0] = params.N;
      population[i]++;
      for(int j = 0; j<params.N; j++)
        population[i][j] = 0.;
    }

    // initialize newed space
    for(int i=0; i<params.lambda; ++i)
    {
      functionValues[i] = std::numeric_limits<T>::max();
    }
    for(int i=0; i<historySize; ++i)
    {
      funcValueHistory[i] = std::numeric_limits<T>::max();
    }
    for(int i=0; i<params.N; ++i)
      for(int j = 0; j<i; ++j)
        C[i][j] = B[i][j] = B[j][i] = 0.;

    for(int i=0; i<params.N; ++i)
    {
      B[i][i] = T(1);
      C[i][i] = rgD[i] = params.rgInitialStds[i]*std::sqrt(params.N/trace);
      C[i][i] *= C[i][i];
      pc[i] = ps[i] = T(0);
    }
    minEW = minElement(rgD, params.N);
    minEW = minEW*minEW;
    maxEW = maxElement(rgD, params.N);
    maxEW = maxEW*maxEW;

    maxdiagC = C[0][0];
    for(int i=1; i < params.N; ++i) if(maxdiagC < C[i][i]) maxdiagC = C[i][i];
    mindiagC = C[0][0];
    for(int i=1; i < params.N; ++i) if(mindiagC > C[i][i]) mindiagC = C[i][i];

    for(int i=0; i<params.N; ++i)
      xmean[i] = xold[i] = params.xstart[i];
    // use in case xstart as typicalX
    if(params.typicalXcase)
      for(int i=0; i<params.N; ++i)
        xmean[i] += sigma*rgD[i]*rand.gauss();

    if(params.resumefile != "")
      resumeDistribution(params.resumefile);

    return publicFitness;
  }

  /**
   * Well, says hello.
   * @return eg. "(5,10)-CMA-ES(mu_eff=3.4), Ver="1.0alpha", dimension=9"
   */
  std::string sayHello()
  {
    std::stringstream stream;
    stream << "(" << params.mu << "," << params.lambda << ")-CMA-ES(mu_eff="
        << std::setprecision(1) << params.mueff << "), Ver=\"" << version
        << "\", dimension=" << params.N << ", diagonalIterations="
        << (long) params.diagonalCov << " (" << getTimeStr() << ")";
    return stream.str();
  }

  /**
   * Allows to restart with saved internal state (distribution) variables (use
   * writeToFile() for saving). Keyword "resume" followed by a filename in
   * initials.par invokes this function during initialization. Searches in
   * filename for the last occurrence of word "resume", followed by a dimension
   * number, and reads the subsequent values for xmean, evolution paths ps and
   * pc, sigma and covariance matrix.  Note that init() needs to be called
   * before calling resume_distribution() explicitely.  In the former all the
   * remaining (strategy-)parameters are set. It can be useful to edit the
   * written parameters, in particular to increase sigma, before resume.
   *
   * Not all internal state parameters are recovered. In particular generation
   * number and xbestever are not restored. For covariance matrices with large
   * condition numbers the writing precision of 6 digits is not sufficient and
   * resume will lead to poor result.
   * @param filename A file, that was written presumably by writeToFile().
   */
  void resumeDistribution(const std::string& filename)
  {
    std::ifstream file(filename);
    if(!file.is_open())
      throw std::runtime_error("resumeDistribution(): could not open '" + filename + "'");

    std::streampos lastResume = 0;
    std::string entry = "";
    while(!file.eof())
    {
      file >> entry;
      if(entry == "resume")
      {
        lastResume = file.tellg();
        break;
      }
    }
    file.clear();
    file.seekg(lastResume);

    int n = 0;
    file >> n;
    if(n != params.N)
      throw std::runtime_error("resumeDistribution(): Dimension numbers do not match");

    // find next "xmean" entry
    while(!file.eof())
    {
      file >> entry;
      if(entry == "xmean")
        break;
    }
    // read xmean
    if(file.eof())
      throw std::runtime_error("resumeDistribution(): 'xmean' not found");
    for(int i=0; i<n; ++i)
      file >> xmean[i];
    file.clear();
    file.seekg(lastResume);

    // find next "path for sigma" entry
    while(!file.eof())
    {
      file >> entry;
      if(entry == "path")
      {
        std::string temp = "";
        file >> temp;
        entry += " " + temp;
        file >> temp;
        entry += " " + temp;
        if(entry == "path for sigma")
          break;
      }
    }
    // read ps
    if(file.eof())
      throw std::runtime_error("resumeDistribution(): 'path for sigma' not found");
    for(int i=0; i<n; ++i)
      file >> ps[i];
    file.clear();
    file.seekg(lastResume);

    // find next "path for C" entry
    while(!file.eof())
    {
      file >> entry;
      if(entry == "path")
      {
        std::string temp = "";
        file >> temp;
        entry += " " + temp;
        file >> temp;
        entry += " " + temp;
        if(entry == "path for C")
          break;
      }
    }
    // read pc
    if(file.eof())
      throw std::runtime_error("resumeDistribution(): 'path for C' not found");
    for(int i=0; i<n; ++i)
      file >> pc[i];
    file.clear();
    file.seekg(lastResume);

    // find next "sigma" entry
    while(!file.eof())
    {
      file >> entry;
      if(entry == "sigma")
        break;
    }
    // read pc
    if(file.eof())
      throw std::runtime_error("resumeDistribution(): 'sigma' not found");
    file >> sigma;
    file.clear();
    file.seekg(lastResume);

    // find next "covariance matrix" entry
    while(!file.eof())
    {
      file >> entry;
      if(entry == "covariance")
      {
        std::string temp = "";
        file >> temp;
        entry += " " + temp;
        if(entry == "covariance matrix")
          break;
      }
    }
    // read C
    if(file.eof())
      throw std::runtime_error("resumeDistribution(): 'covariance matrix' not found");
    for(int i=0; i<params.N; ++i)
      for(int j = 0; j <= i; ++j)
        file >> C[i][j];

    eigensysIsUptodate = false;
    isResumeDone = true;
    updateEigensystem(true);
  }

  /**
   * The search space vectors have a special form: they are arrays with N+1
   * entries. Entry number -1 is the dimension of the search space N.
   * @return A pointer to a "population" of lambda N-dimensional multivariate
   * normally distributed samples.
   */
  T* const* samplePopulation()
  {
    bool diag = params.diagonalCov == 1 || params.diagonalCov >= gen;

    // calculate eigensystem
    if(!eigensysIsUptodate)
    {
      if(!diag)
        updateEigensystem(false);
      else
      {
        for(int i=0; i<params.N; ++i)
          rgD[i] = std::sqrt(C[i][i]);
        minEW = square(minElement(rgD, params.N));
        maxEW = square(maxElement(rgD, params.N));
        eigensysIsUptodate = true;
        eigenTimings.start();
      }
    }

    testMinStdDevs();

    for(int iNk = 0; iNk < params.lambda; ++iNk)
    { // generate scaled random vector D*z
      T* rgrgxink = population[iNk];
      for(int i=0; i<params.N; ++i)
        if(diag)
          rgrgxink[i] = xmean[i] + sigma*rgD[i]*rand.gauss();
        else
          tempRandom[i] = rgD[i]*rand.gauss();
      if(!diag)
        for(int i=0; i<params.N; ++i) // add mutation sigma*B*(D*z)
        {
          T sum = 0.;
          for(int j = 0; j<params.N; ++j)
            sum += B[i][j]*tempRandom[j];
          rgrgxink[i] = xmean[i] + sigma*sum;
        }
    }

    if(state == UPDATED || gen == 0)
      ++gen;
    state = SAMPLED;

    return population;
  }

  /**
   * Can be called after samplePopulation() to resample single solutions of the
   * population as often as desired. Useful to implement a box constraints
   * (boundary) handling.
   * @param i Index to an element of the returned value of samplePopulation().
   *          population[index] will be resampled where \f$0\leq i<\lambda\f$
   *          must hold.
   * @return A pointer to the resampled "population".
   */
  T* const* reSampleSingle(int i)
  {
    T* x;
    assert(i >= 0 && i < params.lambda &&
        "reSampleSingle(): index must be between 0 and sp.lambda");
    x = population[i];
    addMutation(x);
    return population;
  }

  /**
   * Can be called after samplePopulation() to resample single solutions. In
   * general, the function can be used to sample as many independent
   * mean+sigma*Normal(0,C) distributed vectors as desired.
   *
   * Input x can be a pointer to an element of the vector returned by
   * samplePopulation() but this is inconsistent with the const qualifier of the
   * returned value and therefore rather reSampleSingle() should be used.
   * @param x Solution vector that gets sampled a new value. If x == NULL new
   *          memory is allocated and must be released by the user using
   *          delete[].
   * @return A pointer to the resampled solution vector, equals input x for
   *         x != NULL on input.
   */
  T* sampleSingleInto(T* x)
  {
    if(!x)
      x = new T[params.N];
    addMutation(x);
    return x;
  }

  /**
   * Can be called after samplePopulation() to resample single solutions. In
   * general, the function can be used to sample as many independent
   * mean+sigma*Normal(0,C) distributed vectors as desired.
   * @param x Element of the return value of samplePopulation(), that is
   *          pop[0..\f$\lambda\f$]. This solution vector of the population gets
   *          sampled a new value.
   * @return A pointer to the resampled "population" member.
   */
  T const* reSampleSingleOld(T* x)
  {
    assert(x && "reSampleSingleOld(): Missing input x");
    addMutation(x);
    return x;
  }

  /**
   * Used to reevaluate a slightly disturbed solution for an uncertaintly
   * measurement. In case if x == NULL on input, the memory of the returned x
   * must be released.
   * @param x Solution vector that gets sampled a new value. If x == NULL new
   *          memory is allocated and must be released by the user using
   *          delete[] x.
   * @param pxmean Mean vector \f$\mu\f$ for perturbation.
   * @param eps Scale factor \f$\epsilon\f$ for perturbation:
   *            \f$x \sim \mu + \epsilon \sigma N(0,C)\f$.
   * @return A pointer to the perturbed solution vector, equals input x for
   *         x != NULL.
   */
  T* perturbSolutionInto(T* x, T const* pxmean, T eps)
  {
    if(!x)
      x = new T[params.N];
    assert(pxmean && "perturbSolutionInto(): pxmean was not given");
    addMutation(x, eps);
    return x;
  }

  /**
   * Core procedure of the CMA-ES algorithm. Sets a new mean value and estimates
   * the new covariance matrix and a new step size for the normal search
   * distribution.
   * @param fitnessValues An array of \f$\lambda\f$ function values.
   * @return Mean value of the new distribution.
   */
  T* updateDistribution(const T* fitnessValues)
  {
    const int N = params.N;
    bool diag = params.diagonalCov == 1 || params.diagonalCov >= gen;

    assert(state != UPDATED && "updateDistribution(): You need to call "
          "samplePopulation() before update can take place.");
    assert(fitnessValues && "updateDistribution(): No fitness function value array input.");

    if(state == SAMPLED) // function values are delivered here
      countevals += params.lambda;
    else if(params.logWarnings)
      params.logStream <<  "updateDistribution(): unexpected state" << std::endl;

    // assign function values
    for(int i=0; i<params.lambda; ++i)
      population[i][N] = functionValues[i] = fitnessValues[i];

    // Generate index
    sortIndex(fitnessValues, index, params.lambda);

    // Test if function values are identical, escape flat fitness
    if(fitnessValues[index[0]] == fitnessValues[index[(int) params.lambda / 2]])
    {
      sigma *= std::exp(T(0.2) + params.cs / params.damps);
      if(params.logWarnings)
      {
        params.logStream << "Warning: sigma increased due to equal function values"
            << std::endl << "   Reconsider the formulation of the objective function";
      }
    }

    // update function value history
    for(int i = (int) *(funcValueHistory - 1) - 1; i > 0; --i)
      funcValueHistory[i] = funcValueHistory[i - 1];
    funcValueHistory[0] = fitnessValues[index[0]];

    // update xbestever
    if(xBestEver[N] > population[index[0]][N] || gen == 1)
      for(int i=0; i <= N; ++i)
      {
        xBestEver[i] = population[index[0]][i];
        xBestEver[N+1] = countevals;
      }

    const T sqrtmueffdivsigma = std::sqrt(params.mueff) / sigma;
    // calculate xmean and rgBDz~N(0,C)
    for(int i=0; i<N; ++i)
    {
      xold[i] = xmean[i];
      xmean[i] = 0.;
      for(int iNk = 0; iNk < params.mu; ++iNk)
        xmean[i] += params.weights[iNk]*population[index[iNk]][i];
      BDz[i] = sqrtmueffdivsigma*(xmean[i]-xold[i]);
    }

    // calculate z := D^(-1)* B^(-1)* rgBDz into rgdTmp
    for(int i=0; i<N; ++i)
    {
      T sum;
      if(diag)
        sum = BDz[i];
      else
      {
        sum = 0.;
        for(int j = 0; j<N; ++j)
          sum += B[j][i]*BDz[j];
      }
      tempRandom[i] = sum/rgD[i];
    }

    // cumulation for sigma (ps) using B*z
    const T sqrtFactor = std::sqrt(params.cs*(T(2)-params.cs));
    const T invps = T(1)-params.cs;
    for(int i=0; i<N; ++i)
    {
      T sum;
      if(diag)
        sum = tempRandom[i];
      else
      {
        sum = T(0);
        T* Bi = B[i];
        for(int j = 0; j<N; ++j)
          sum += Bi[j]*tempRandom[j];
      }
      ps[i] = invps*ps[i] + sqrtFactor*sum;
    }

    // calculate norm(ps)^2
    T psxps(0);
    for(int i=0; i<N; ++i)
    {
      const T& rgpsi = ps[i];
      psxps += rgpsi*rgpsi;
    }

    // cumulation for covariance matrix (pc) using B*D*z~N(0,C)
    int hsig = std::sqrt(psxps) / std::sqrt(T(1) - std::pow(T(1) - params.cs, T(2)* gen))
        / chiN < T(1.4) + T(2) / (N + 1);
    const T ccumcovinv = 1.-params.ccumcov;
    const T hsigFactor = hsig*std::sqrt(params.ccumcov*(T(2)-params.ccumcov));
    for(int i=0; i<N; ++i)
      pc[i] = ccumcovinv*pc[i] + hsigFactor*BDz[i];

    // update of C
    adaptC2(hsig);

    // update of sigma
    sigma *= std::exp(((std::sqrt(psxps) / chiN) - T(1))* params.cs / params.damps);

    state = UPDATED;
    return xmean;
  }

  /**
   * Request a scalar parameter from CMA-ES.
   * @param key Key of the requested scalar.
   * @return The desired value.
   */
  T get(GetScalar key)
  {
    switch(key)
    {
      case AxisRatio:     return maxElement(rgD, params.N) / minElement(rgD, params.N);
      case Eval:          return countevals;
      case Fitness:       return functionValues[index[0]];
      case FBestEver:     return xBestEver[params.N];
      case Generation:    return gen;
      case MaxEval:       return params.stopMaxFunEvals;
      case MaxIter:       return std::ceil(params.stopMaxIter);
      case MaxAxisLength: return sigma*std::sqrt(maxEW);
      case MinAxisLength: return sigma*std::sqrt(minEW);
      case MaxStdDev:     return sigma*std::sqrt(maxdiagC);
      case MinStdDev:     return sigma*std::sqrt(mindiagC);
      case Dimension:     return params.N;
      case SampleSize:    return params.lambda;
      case Sigma:         return sigma;
      default:       throw std::runtime_error("get(): No match found for key");
    }
  }

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @return Pointer to the desired value array. Its content might be
   *         overwritten during the next call to any member functions other
   *         than get().
   */
  const T* getPtr(GetVector key)
  {
    switch(key)
    {
      case DiagC:
      {
        for(int i=0; i<params.N; ++i)
          output[i] = C[i][i];
        return output;
      }
      case DiagD:
        return rgD;
      case StdDev:
      {
        for(int i=0; i<params.N; ++i)
          output[i] = sigma*std::sqrt(C[i][i]);
        return output;
      }
      case XBestEver:
        return xBestEver;
      case XBest:
        return population[index[0]];
      case XMean:
        return xmean;
      default:
        throw std::runtime_error("getPtr(): No match found for key");
    }
  }

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @return Pointer to the desired value array with unlimited reading and
   *         writing access to its elements. The memory must be explicitly
   *         released using delete[].
   */
  T* getNew(GetVector key)
  {
    return getInto(key, 0);
  }

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @param res Memory of size N == dimension, where the desired values are
   *            written into. For mem == NULL new memory is allocated as with
   *            calling getNew() and must be released by the user at some point.
   */
  T* getInto(GetVector key, T* res)
  {
    T const* res0 = getPtr(key);
    if(!res)
      res = new T[params.N];
    for(int i=0; i<params.N; ++i)
      res[i] = res0[i];
    return res;
  }

  /**
   * Some stopping criteria can be set in initials.par, with names starting
   * with stop... Internal stopping criteria include a maximal condition number
   * of about 10^15 for the covariance matrix and situations where the numerical
   * discretisation error in x-space becomes noticeably. You can get a message
   * that contains the matched stop criteria via getStopMessage().
   * @return Does any stop criterion match?
   */
  bool testForTermination()
  {
    T range, fac;
    int iAchse, iKoo;
    int diag = params.diagonalCov == 1 || params.diagonalCov >= gen;
    int N = params.N;
    std::stringstream message;

    if(stopMessage != "")
    {
      message << stopMessage << std::endl;
    }

    // function value reached
    if((gen > 1 || state > SAMPLED) && params.stStopFitness.flg &&
        functionValues[index[0]] <= params.stStopFitness.val)
    {
      message << "Fitness: function value " << functionValues[index[0]]
          << " <= stopFitness (" << params.stStopFitness.val << ")" << std::endl;
    }

    // TolFun
    range = std::max(maxElement(funcValueHistory, (int) std::min(gen, *(funcValueHistory - 1))),
        maxElement(functionValues, params.lambda)) -
        std::min(minElement(funcValueHistory, (int) std::min(gen, *(funcValueHistory - 1))),
        minElement(functionValues, params.lambda));

    if(gen > 0 && range <= params.stopTolFun)
    {
      message << "TolFun: function value differences " << range
          << " < stopTolFun=" << params.stopTolFun << std::endl;
    }

    // TolFunHist
    if(gen > *(funcValueHistory - 1))
    {
      range = maxElement(funcValueHistory, (int) *(funcValueHistory - 1))
          - minElement(funcValueHistory, (int) *(funcValueHistory - 1));
      if(range <= params.stopTolFunHist)
        message << "TolFunHist: history of function value changes " << range
            << " stopTolFunHist=" << params.stopTolFunHist << std::endl;
    }

    // TolX
    int cTemp = 0;
    for(int i=0; i<N; ++i)
    {
      cTemp += (sigma*std::sqrt(C[i][i]) < params.stopTolX) ? 1 : 0;
      cTemp += (sigma*pc[i] < params.stopTolX) ? 1 : 0;
    }
    if(cTemp == 2*N)
    {
      message << "TolX: object variable changes below " << params.stopTolX << std::endl;
    }

    // TolUpX
    for(int i=0; i<N; ++i)
    {
      if(sigma*std::sqrt(C[i][i]) > params.stopTolUpXFactor*params.rgInitialStds[i])
      {
        message << "TolUpX: standard deviation increased by more than "
            << params.stopTolUpXFactor << ", larger initial standard deviation recommended."
            << std::endl;
        break;
      }
    }

    // Condition of C greater than dMaxSignifKond
    if(maxEW >= minEW* dMaxSignifKond)
    {
      message << "ConditionNumber: maximal condition number " << dMaxSignifKond
          << " reached. maxEW=" << maxEW <<  ",minEW=" << minEW << ",maxdiagC="
          << maxdiagC << ",mindiagC=" << mindiagC << std::endl;
    }

    // Principal axis i has no effect on xmean, ie. x == x + 0.1* sigma* rgD[i]* B[i]
    if(!diag)
    {
      for(iAchse = 0; iAchse < N; ++iAchse)
      {
        fac = 0.1* sigma* rgD[iAchse];
        for(iKoo = 0; iKoo < N; ++iKoo)
        {
          if(xmean[iKoo] != xmean[iKoo] + fac* B[iKoo][iAchse])
            break;
        }
        if(iKoo == N)
        {
          message << "NoEffectAxis: standard deviation 0.1*" << (fac / 0.1)
              << " in principal axis " << iAchse << " without effect" << std::endl;
          break;
        }
      }
    }
    // Component of xmean is not changed anymore
    for(iKoo = 0; iKoo < N; ++iKoo)
    {
      if(xmean[iKoo] == xmean[iKoo] + sigma*std::sqrt(C[iKoo][iKoo])/T(5))
      {
        message << "NoEffectCoordinate: standard deviation 0.2*"
            << (sigma*std::sqrt(C[iKoo][iKoo])) << " in coordinate " << iKoo
            << " without effect" << std::endl;
        break;
      }
    }

    if(countevals >= params.stopMaxFunEvals)
    {
      message << "MaxFunEvals: conducted function evaluations " << countevals
          << " >= " << params.stopMaxFunEvals << std::endl;
    }
    if(gen >= params.stopMaxIter)
    {
      message << "MaxIter: number of iterations " << gen << " >= "
          << params.stopMaxIter << std::endl;
    }

    stopMessage = message.str();
    return stopMessage != "";
  }

  /**
   * A message that contains a detailed description of the matched stop
   * criteria.
   */
  std::string getStopMessage()
  {
    return stopMessage;
  }

  /**
   * @param filename Output file name.
   * @param key Key of type WriteKey that indicates the content that should be
   *            written. You can combine multiple keys with |.
   */
  void writeToFile(int key, const std::string& filename)
  {
    std::ofstream file(filename, std::ios_base::app);

    if(file.is_open())
    {
      if(gen > 0 || filename.substr(0, 11) != "outcmaesfit")
        writeToStream(key, file); /* do not write fitness for gen==0 */
      file.close();
    }
    else
    {
      throw std::runtime_error("writeToFile(): could not open '" + filename + "'");
    }
  }

  /**
   * Conducts the eigendecomposition of C into B and D such that
   * \f$C = B \cdot D \cdot D \cdot B^T\f$ and \f$B \cdot B^T = I\f$
   * and D diagonal and positive.
   * @param force For force == true the eigendecomposion is conducted even if
   *              eigenvector and values seem to be up to date.
   */
  void updateEigensystem(bool force)
  {
    eigenTimings.update();

    if(!force)
    {
      if(eigensysIsUptodate)
        return;
      // return on modulo generation number
      if(gen < genOfEigensysUpdate + params.updateCmode.modulo)
        return;
      // return on time percentage
      if(params.updateCmode.maxtime < 1.00
          && eigenTimings.tictoctime > params.updateCmode.maxtime* eigenTimings.totaltime
          && eigenTimings.tictoctime > 0.0002)
        return;
    }

    eigenTimings.tic();
    eigen(rgD, B, tempRandom);
    eigenTimings.toc();

    // find largest and smallest eigenvalue, they are supposed to be sorted anyway
    minEW = minElement(rgD, params.N);
    maxEW = maxElement(rgD, params.N);

    if(doCheckEigen) // needs O(n^3)! writes, in case, error message in error file
      checkEigen(rgD, B);

    for(int i=0; i<params.N; ++i)
      rgD[i] = std::sqrt(rgD[i]);

    eigensysIsUptodate = true;
    genOfEigensysUpdate = gen;
  }

  /**
   * Distribution mean could be changed before samplePopulation(). This might
   * lead to unexpected behaviour if done repeatedly.
   * @param newxmean new mean, if it is NULL, it will be set to the current mean
   * @return new mean
   */
  T const* setMean(const T* newxmean)
  {
    assert(state != SAMPLED && "setMean: mean cannot be set inbetween the calls"
        "of samplePopulation and updateDistribution");

    if(newxmean && newxmean != xmean)
      for(int i=0; i<params.N; ++i)
        xmean[i] = newxmean[i];
    else
      newxmean = xmean;

    return newxmean;
  }
};
