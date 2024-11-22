#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace arma ;

// [[Rcpp::export]]
arma::mat get_sandwich_variance_B(
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & A,
  const arma::mat & S,
  const arma::mat & Sigma,
  const arma::vec & Diag_Omega
) {

  const uword n = Y.n_rows ;
  const uword p = Y.n_cols ;
  const uword d = X.n_cols ;

  mat Cn(d*p, d*p, fill::zeros) ;
  mat Dn(d*p, d*p, fill::zeros) ;

  mat D = pow(A, -1) + pow(S, 4) / (1 + square(S) % (A + ones(n) * Diag_Omega.t())) ;
  mat YmA = Y - A ;

  for (uword i=0; i<n; i++) {
    mat xxt_i = X.row(i).t() * X.row(i) ;
    mat yyt_i = YmA.row(i).t() * YmA.row(i) ;
    Cn -= kron(inv_sympd(Sigma + diagmat(D.row(i))), xxt_i) / n ;
    Dn += kron(yyt_i, xxt_i) / n ;
  }

  mat Cn_inv = inv(Cn) ;

  return (Cn_inv * Dn * Cn_inv) / n ;
}
