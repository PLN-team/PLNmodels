#ifndef _data_struct_H
#define _data_struct_H

#include <RcppArmadillo.h>

inline arma::mat logfact(arma::mat Y) {
  arma::mat v = Y.replace(0, 1);
  return sum(v % arma::log(v) - v + arma::log(8*pow(v,3) + 4*pow(v, 2) + v + 1/30)/6 + std::log(M_PI)/2, 1);
}

typedef std::vector<double> stdvec;

typedef struct optim_data {
    arma::mat Y          ;
    arma::mat X          ;
    arma::mat O          ;
    arma::vec w          ;
    arma::mat Omega      ;
    arma::mat Theta      ;
    double w_bar         ;
    double log_det_Omega ;
    arma::vec Ki         ;
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
               const arma::mat &offsets,
               const arma::mat &weights
    ) : Y(responses), X(covariates), O(offsets), w(weights)
      {
        n = Y.n_rows ;
        p = Y.n_cols ;
        d = X.n_cols ;
        iterations = 0     ;
        w_bar = accu(w)    ;
        Ki = - logfact(Y) + .5 * (1+(1-p)* std::log(2*M_PI)) ;
      } ;

} optim_data ;

#endif
