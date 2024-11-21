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


  arma::mat get_iCnB = [&A, &S, &D_omega, &Sigma](
      ) {


    return ;
  }



  return ;
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
