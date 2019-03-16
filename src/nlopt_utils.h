#ifndef _nlopt_utils_H
#define _nlopt_utils_H

#include "RcppArmadillo.h"
#include "data_struct.h"
#include "nlopt.h"

// Convert string to nlopt_alogirthm
nlopt_algorithm getAlgorithmCode(const std::string &) ;

// Initialize the nlopt optimizer
nlopt_opt initNLOPT(int, Rcpp::List) ;

#endif
