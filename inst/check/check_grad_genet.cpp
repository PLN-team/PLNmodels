#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

// [[Rcpp::export]]
double objective_genet(arma::vec param, const arma::mat &Y, const arma::mat&X, const arma::mat&V, const arma::vec &Lambda) {
  arma::uword p = Y.n_cols;
  arma::uword d = X.n_cols;
  arma::uword n = Y.n_rows;

  arma::mat Theta(&param[0], p, d, true);
  arma::mat M(&param[p*d], n, p, true);
  arma::mat S(&param[p*d + n*p], n, p, true);
  double rho  = param[p*d + 2*n*p] ;

  arma::mat S2 = S % S;
  arma::mat Z = X * Theta.t() + M;
  arma::mat A = exp(Z + 0.5 * S2);
  arma::vec u = rho * Lambda + (1-rho) ;
  arma::mat R = V.t() * (M.t() * M + diagmat(arma::ones(n).t() * S2)) * V ;
  double sigma2 = accu( diagvec(R) / u ) / (double(p * n));
  arma::mat Omega = V * diagmat(pow(sigma2 * u, -1)) * V.t() ;

  double objective = accu(A - Y % Z - 0.5 * log(S2)) +
    //          0.5 * n * accu(log(u * sigma2));
     0.5 * trace(Omega * (M.t() * M + diagmat(arma::ones(n).t() * S2))) +
     0.5 * n * accu(log(u * sigma2));

  return objective ;
};

// [[Rcpp::export]]
arma::vec grad_genet(arma::vec param, const arma::mat &Y, const arma::mat&X, const arma::mat&V, const arma::vec &Lambda) {
  arma::uword p = Y.n_cols;
  arma::uword d = X.n_cols;
  arma::uword n = Y.n_rows;

  arma::mat Theta(&param[0], p, d, true);
  arma::mat M(&param[p*d], n, p, true);
  arma::mat S(&param[p*d + n*p], n, p, true);
  double rho  = param[p*d + 2*n*p] ;

  arma::mat S2 = S % S;
  arma::mat Z = X * Theta.t() + M;
  arma::mat A = exp(Z + 0.5 * S2);
  arma::vec u = rho * Lambda + (1-rho) ;
  arma::mat R = V.t() * (M.t() * M + diagmat(arma::ones(n).t() * S2)) * V ;
  double sigma2 = accu( diagvec(R) / u ) / (double(p * n));
  arma::mat Omega = V * diagmat(pow(sigma2 * u, -1)) * V.t() ;

  arma::vec Theta_grad = arma::vectorise((A - Y).t() * X);
  arma::vec M_grad =  arma::vectorise(M * Omega + A - Y);
  arma::vec S_grad =  arma::vectorise(S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));
  arma::vec rho_grad(1);
  rho_grad[0] = accu(0.5 * n * (Lambda - 1) / u - (0.5/sigma2) * diagvec(R) % (Lambda - 1) / pow(u, 2) ) ;

  arma::vec grad = join_cols(Theta_grad, M_grad, S_grad, rho_grad);

  return grad ;
};


/*** R
library(nloptr)
library(mvtnorm)
n <- 100; p <- 5; d <- 1
corr <- toeplitz(0.5^(1:p - 1))
sigma2 <- 2; rho <- 0.25
Sigma <- sigma2 * (rho * corr + (1 - rho) * diag(1, p, p))
Z <- rmvnorm(n, sigma = corr)
Y <- matrix(rpois(n* p, exp(Z)), n, p)
X <- matrix(1,n,d)

eig <- eigen(corr)
param <- c(rep(0, p*d) , c(matrix(0, n, p)), c(matrix(sqrt(0.1), n, p)), 0.5)
print(objective_genet(param, Y, X, eig$vectors, eig$values))
print(grad_genet(param, Y, X, eig$vectors, eig$values))

objective_genet_R<- function(.x) {
  objective_genet(.x, Y = Y, X = X, V = eig$vectors, Lambda = eig$values)
}

grad_genet_R<- function(.x) {
  grad_genet(.x, Y = Y, X = X, V = eig$vectors, Lambda = eig$values)
}
check <- nloptr::check.derivatives(.x = param, objective_genet_R, grad_genet_R)
*/
