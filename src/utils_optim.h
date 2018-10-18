#ifndef _utils_optim_H
#define _utils_optim_H

#include "RcppArmadillo.h"
#include <nlopt.hpp>

arma::vec get_KY(arma::mat Y, arma::vec w) ;

arma::mat logfact(arma::mat Y) ;

typedef std::vector<double> stdvec;

typedef struct optim_data {
    arma::mat Y          ;
    arma::mat X          ;
    arma::mat O          ;
    arma::vec w          ;
    arma::mat Omega      ;
    arma::mat Theta      ;
    double log_det_Omega ;
    double KY            ;
    arma::vec KYi        ;
    int iterations       ;
    int n                ;
    int p                ;
    int d                ;
    int q                ;

    // constructors

    // Empty constructor
    optim_data() { } ;

    // PLN constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets
    ) : Y(responses), X(covariates), O(offsets)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        iterations = 0 ;
        KYi = logfact(Y);
        KY  = accu(KYi) ;
      } ;
    // weighted PLN constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets,
               const arma::mat &weights
    ) : Y(responses), X(covariates), O(offsets), w(weights)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        iterations = 0 ;
        KYi = logfact(Y) ;
        KY = accu(w % KYi) ;
      } ;
    // PLNPCA constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets,
               const int rank
    ) : Y(responses), X(covariates), O(offsets), q(rank)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        iterations = 0 ;
        KYi = logfact(Y) ;
        KY = accu(KYi) ;
      } ;
    // PLNnetwork constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets,
               const arma::mat covinv,
               const double log_det
    ) : Y(responses), X(covariates), O(offsets), Omega(covinv), log_det_Omega(log_det)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        iterations = 0 ;
        KYi = logfact(Y) ;
        KY = accu(KYi) ;
      } ;
    // PLN VE-step constructor
    optim_data(const arma::mat &responses,
               const arma::mat &covariates,
               const arma::mat &offsets,
               const arma::mat &regression_parameters,
               const arma::mat &covinv,
               const double log_det
    ) : Y(responses), X(covariates), O(offsets), Omega(covinv), Theta(regression_parameters), log_det_Omega(log_det)
    {
      n = Y.n_rows ;
      p = Y.n_cols ;
      d = X.n_cols ;
      iterations = 0 ;
      KYi = logfact(Y) ;
      KY = accu(KYi) ;
    } ;

} optim_data ;

// Convert string to nlopt_alogirthm
nlopt::algorithm getAlgorithmCode( const std::string) ;

nlopt::opt initNLOPT(int, Rcpp::List) ;

#endif
