#ifndef _utils_optim_H
#define _utils_optim_H

#include "RcppArmadillo.h"
#include <nlopt.hpp>

double K(arma::mat Y) ;

typedef std::vector<double> stdvec;

typedef struct optim_data {
    arma::mat Y    ;
    arma::mat X    ;
    arma::mat O    ;
    double KY      ;
    int iterations ;
    int n          ;
    int p          ;
    int d          ;
    int q          ;

    // constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets,
               const int rank= -1) :
      Y(responses), X(covariates), O(offsets), q(rank)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        if (rank == -1) {
          q = p ;
        } else {
          q = rank ;
        }
        iterations = 0 ;
        KY = K(Y) ;

      } ;
} optim_data ;

// Convert string to nlopt_alogirthm
nlopt::algorithm getAlgorithmCode( const std::string) ;

nlopt::opt initNLOPT(int, Rcpp::List) ;

#endif
