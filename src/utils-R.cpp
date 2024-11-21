#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

// [[Rcpp::export]]
arma::mat get_sandwich_variance_B(
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & A,
  const arma::mat & S,
  const arma::mat & Sigma,
  const arma::vec & Diag_Omega
) {

  arma::uword n = Y.n_rows ;
  arma::uword p = Y.n_cols ;
  arma::uword d = X.n_cols ;

  auto get_iCnB = [&A, &S, &Diag_Omega, &Sigma](arma::uword i) {
    arma::vec a = A.row(i) ;
    arma::vec s = S.row(i) ;
    arma::mat D = diagmat(pow(a, -1) + pow(s, 4) / (1 + pow(s,2) * (a + Diag_Omega))) ;

    return arma::inv_sympd(Sigma + D);
  };

  arma::mat YmA = Y - A ;
  arma::mat Cn = arma::zeros(d*p, d*p) ;
  arma::mat Dn = arma::zeros(d*p, d*p) ;
  for (int i=0; i<n; i++) {
    arma::mat xxt_i = X.col(i) * X.col(i).t() ;
    arma::mat yyt_i = YmA.col(i) * YmA.col(i).t() ;
    Cn = Cn - arma::kron(get_iCnB(i), xxt_i) / n ;
    Dn = Dn + arma::kron(yyt_i, xxt_i) / n ;
  }

  arma::mat Cn_inv = arma::inv_sympd(Cn) ;

  return (Cn_inv * Dn * Cn_inv) / n ;
}

// vcov_sandwich_B = function(Y, X) {
//  getMat_iCnB <- function(i) {
//     a_i <- as.numeric(private$A[i, ])
//     s_i <- as.numeric(private$S[i, ])
//     omega <- as.numeric(diag(private$Omega))
//     diag_mat_i <- diag(1/a_i + s_i^4 / (1 + s_i^2 * (a_i + omega)))
//     solve(private$Sigma + diag_mat_i)
//  }
//
//  YmA <- Y - private$A
//    Dn <- matrix(0, self$d*self$p, self$d*self$p)
//    Cn <- matrix(0, self$d*self$p, self$d*self$p)
//    for (i in 1:self$n) {
//    xxt_i <- tcrossprod(X[i, ])
//    Cn <- Cn - kronecker(getMat_iCnB(i) , xxt_i) / self$n
//    Dn <- Dn + kronecker(tcrossprod(YmA[i,]), xxt_i) / self$n
//  }
//   Cn_inv <- solve(Cn)
//     dim_names <- dimnames(attr(private$B, "vcov_variational"))
//     vcov_sand <- ((Cn_inv %*% Dn %*% Cn_inv) / self$n) %>% `dimnames<-`(dim_names)
//
