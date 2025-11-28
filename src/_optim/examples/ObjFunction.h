/*
 * ObjFunction.hpp
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 */

#ifndef OBJFUNCTION_HPP_
#define OBJFUNCTION_HPP_
#include <algorithm>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>


// Computes the Himmelblau function.
// This function has 4 global minima:
//   X* = (  3,        2       ), F(X*) = 0.
//   X* = (  3.58439, -1.84813 ), F(X*) = 0.
//   X* = ( -3.77934, -3.28317 ), F(X*) = 0.
//   X* = ( -2.80512,  3.13134 ), F(X*) = 0.
// Suggested starting xx are   (+1,+1),   (-1,+1),   (-1,-1),   (+1,-1),
double himmelblau(double *xx, int nx) {
    double sum = 0;
    sum += std::pow(std::pow(xx[0], 2) + xx[1] - 11, 2);
    sum += std::pow(xx[0] + std::pow(xx[1], 2) - 7, 2);
    return sum;
}

// Computes the Rosenbrock function.
// There is a global minimum at X* = (1,1), F(X*) = 0.
// The starting point X = [ -1.2, 1. ] is suggested.
// The contours are sharply twisted.
double rosenbrock(double *xx, int nx) {

    double sum = 0;
    sum += std::pow(1. - xx[0], 2);
    sum += 100. * std::pow(xx[1] - xx[0] * xx[0], 2);

    return sum;
}

// Computes the Beale function.
// This function has a global minimizer:
//   X = ( 3., 0.5 ) for which   F(X) = 0.
// For an easy computation, start the search at   X = ( 1., 1. )
// A harder computation starts at                 X = ( 1., 4. )
double beale(double *xx, int nx) {

    double sum = 0;
    sum += std::pow(1.5   - xx[0] * (1. - xx[1]), 2);
    sum += std::pow(2.25  - xx[0] * (1. - std::pow(xx[1], 2)), 2);
    sum += std::pow(2.625 - xx[0] * (1. - std::pow(xx[1], 3)), 2);

    return sum;

}

// Evaluates the Bohachevsky function #2.
// The minimizer is   X* = [ 0., 0. ]   F(X*) = 0.
// Suggested starting point:   X init = [ 0.6, 1.3 ];
double bohach2(double *xx, int nx) {

    double sum = 0;
    const double pi = 3.14159265358979323846;

    sum += xx[0] * xx[0];
    sum +=  2. * xx[1] * xx[1];
    sum -=  0.3 * cos(3. * pi * xx[0]) * cos(4. * pi * xx[1]);
    sum += + 0.3;

    return sum;

}

// Extended Rosenbrock function.
// The number of dimensions is arbitrary, except that it must be even.
// There is a global minimum at X* = (1,1,...), F(X*) = 0.
// The contours are sharply twisted.
double extended_rosenbrock(double *xx, int nx) {
    if(nx % 2) { std::cerr << "Dimension must be an even number for extended Rosenbrock function.";     exit(1); }

    double sum = 0;
    double *r = new double[nx];

    for (int i = 0; i < nx; i += 2) {
        r[i] = 1. - xx[i];
        r[i + 1] = 10. * (xx[i + 1] - std::pow(xx[i], 2));
    }

   for (int i = 0; i < nx; i++)    sum += std::pow(r[i], 2);

    delete[] r;
    return sum;
}

// Computes the Powell singular quartic function.
// This function has a global minimizer:   X* = ( 0., 0., 0., 0. ), F(X*) = 0.
// Start the search at   X = ( 3., -1., 0., 1. )
double powell(double *xx, int nx) {
    double f1 = xx[0] + 10. * xx[1];
    double f2 = xx[2] - xx[3];
    double f3 = xx[1] - 2. * xx[2];
    double f4 = xx[0] - xx[3];
    return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4;
}

//  Evaluates the Goldstein-Price polynomial.
//  The minimizer is  X* = [ 0., -1. ]  F(X*) = 3.
//  Suggested starting point:  X init = [ -0.5, 0.25 ] (easy convergence),   X init = [ -4., 5.00 ] (harder convergence)
double goldstein_price(double *xx, int nx) {
    double a = xx[0] + xx[1] + 1.;
    double b = 19. - 14. * xx[0] + 3. * xx[0] * xx[0] - 14. * xx[1] + 6. * xx[0] * xx[1] + 3. * xx[1] * xx[1];
    double c = 2. * xx[0] - 3. * xx[1];
    double d = 18. - 32. * xx[0] + 12. * xx[0] * xx[0] + 48. * xx[1] - 36. * xx[0] * xx[1] + 27. * xx[1] * xx[1];
    double f = ( 1. + a * a * b ) * ( 30. + c * c * d );
    return f;
}


//  Computes the local function.
//  This function has a local minimizer:    X* = ( 0.2858054412..., 0.2793263206...), F(X*) = 5.9225...
//  and a global minimizer:    X* = ( -21.02667179..., -36.75997872...), F(X*) = 0.
//  Suggested starting point for local minimizer:    X = ( 1, 1 ), F(X) = 3.33 * 10^6.
//  Suggested starting point for global minimizer:   X = ( -15, -35), F(X) = 1.49 * 10^8.
double local(double *xx, int nx) {
    double sum = 0;
    sum += std::pow(std::pow(xx[0], 2) + 12 * xx[1] - 1, 2);
    sum += std::pow(49 * std::pow(xx[0], 2) + 49 * std::pow(xx[1], 2) + 84 * xx[0] + 2324 * xx[1] - 681, 2);
    return sum;
}

//  Computes the sum of squares function
double objFunction1(double *xx, int nx) {
    double sum = 0;
    for (int i = 0; i < nx; ++i) sum += std::pow(xx[i] - 2., 2) / nx;
    return sum;
}

//  Computes the sum of absolute values function
double objFunction2(double *xx, int nx) {
    double sum = 0;
    for (int i = 0; i < nx; ++i)   sum += std::abs(xx[i] - 2.)/ nx;
    return sum;
}

//  Computes the sum of squares function
double objFunction3(double *xx, int nx) {
    double sum = 0;
    for (int i=0; i<nx; ++i)  sum += std::pow(xx[i] - (double)i/2., 2) / nx;
    return sum;
}
//  Computes the sum of squares function
double objFunctionWebPy(double *xx, int nx) {
    //double sum = 0;
    //for (int i=0; i<nx; ++i)
      //sum += .5*pow(1 - xx[i-i%2],2) + pow(xx[i] - pow(xx[i-i%2],2),2) / nx;
    //return sum;

  int i;
  double sum = 1e4*xx[0]*xx[0] + 1e-4*xx[1]*xx[1];
  for(i = 2; i < nx; ++i)
    sum += xx[i]*xx[i];
  return sum;
}

#endif
