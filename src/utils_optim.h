#ifndef _utils_optim_H
#define _utils_optim_H

#include "RcppArmadillo.h"
#include <nlopt.hpp>

typedef std::vector<double> stdvec;

struct optim_data {
    arma::mat Y    ;
    arma::mat X    ;
    arma::mat O    ;
    double KY      ;
    int iterations ;
};

double K(arma::mat Y) ;

// Convert string to nlopt_alogirthm
nlopt::algorithm getAlgorithmCode( const std::string algorithm_str) ;

nlopt::opt initNLOPT(int n_param, Rcpp::List control) ;

#endif
