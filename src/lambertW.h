#include <Rcpp.h>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace Rcpp;

const double EPS = 2.2204460492503131e-16;
const double M_1_E = 1.0 / M_E;

double FritschIter(double x, double w) ;
double lambertW0_CS(double x) ;
double lambertWm1_CS(double x) ;


