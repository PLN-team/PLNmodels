#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
double objective_Theta0(arma::vec param, const arma::mat &X, const arma::mat&Pi) {

  arma::uword n = Pi.n_rows;
  arma::uword p = Pi.n_cols;
  arma::uword d = X.n_cols;

  arma::mat Theta0(&param[0], d, p, true);
  arma::mat Xt_Pi = X.t() * Pi;
  arma::mat e_X_Theta0 = exp(X * Theta0);
  double objective = -trace(Xt_Pi.t() * Theta0) + accu(log(1. + e_X_Theta0));

  return objective ;
}

// [[Rcpp::export]]
arma::vec grad_Theta0(arma::vec param, const arma::mat&X, const arma::mat&Pi) {

  arma::uword n = Pi.n_rows;
  arma::uword p = Pi.n_cols;
  arma::uword d = X.n_cols;

  arma::mat Theta0(&param[0], d, p, true);
  arma::mat Xt_Pi = X.t() * Pi;
  arma::mat e_X_Theta0 = exp(X * Theta0);
  arma::vec grad = arma::vectorise(-Xt_Pi + X.t() * pow(1. + e_X_Theta0, -1));

  return grad ;
};

// _____________________________________________________________________________

// [[Rcpp::export]]
double objective_M(arma::vec param,
                    const arma::mat &Y,
                    const arma::mat &X,
                    const arma::mat &O,
                    const arma::mat &Pi,
                    const arma::mat &S,
                    const arma::mat &Theta,
                    const arma::mat &Omega) {

  arma::uword n = Y.n_rows;
  arma::uword p = Y.n_cols;
  arma::mat M(&param[0], n, p, true);

  const arma::mat X_Theta = X * Theta; // (n,p)
  const arma::mat O_S2 = O + 0.5 * S % S; // (n,p)
  arma::mat A = exp(O_S2 + M);                       // (n,p)
  arma::mat M_X_Theta_Omega = (M - X_Theta) * Omega; // (n,p)

  double objective = trace((Pi - 1.).t() * (Y % M - A)) + 0.5 * trace(M_X_Theta_Omega * (M - X_Theta).t());
  return objective;

}

// [[Rcpp::export]]
arma::vec grad_M(arma::vec param,
                    const arma::mat &Y,
                    const arma::mat &X,
                    const arma::mat &O,
                    const arma::mat &Pi,
                    const arma::mat &S,
                    const arma::mat &Theta,
                    const arma::mat &Omega) {

  arma::uword n = Y.n_rows;
  arma::uword p = Y.n_cols;
  arma::mat M(&param[0], n, p, true);
  const arma::mat X_Theta = X * Theta; // (n,p)
  const arma::mat O_S2 = O + 0.5 * S % S; // (n,p)
  arma::mat A = exp(O_S2 + M);                       // (n,p)
  arma::mat M_X_Theta_Omega = (M - X_Theta) * Omega; // (n,p)

  arma::vec grad = arma::vectorise((Pi - 1.) % (Y - A) + M_X_Theta_Omega);

  return grad ;
};

// [[Rcpp::export]]
double objective_S(arma::vec param,
                    const arma::mat &O,
                    const arma::mat &M,
                    const arma::mat &Pi,
                    const arma::mat &Theta,
                    const arma::mat &Omega) {

  arma::uword n = O.n_rows;
  arma::uword p = O.n_cols;
  arma::mat S(&param[0], n, p, true);

  const arma::mat O_M = O + M;
  const arma::vec diag_Omega = diagvec(Omega);
  arma::mat A = exp(O_M + 0.5 * S % S); // (n,p)

  // trace(1^T log(S)) == accu(log(S)).
  // S_bar = diag(sum(S, 0)). trace(Omega * S_bar) = dot(diagvec(Omega), sum(S2, 0))
  double objective = trace((1. - Pi).t() * A) + 0.5 * dot(diag_Omega, sum(S % S, 0)) - 0.5 * accu(log(S % S));

  return objective;

}

// [[Rcpp::export]]
arma::vec grad_S(arma::vec param,
                    const arma::mat &O,
                    const arma::mat &M,
                    const arma::mat &Pi,
                    const arma::mat &Theta,
                    const arma::mat &Omega) {

  arma::uword n = O.n_rows;
  arma::uword p = O.n_cols;
  arma::mat S(&param[0], n, p, true);

  const arma::mat O_M = O + M;
  const arma::vec diag_Omega = diagvec(Omega);
  arma::mat A = exp(O_M + 0.5 * S % S); // (n,p)

  arma::vec grad = arma::vectorise(S.each_row() % diagvec(Omega).t() + (1. - Pi) % S % A - pow(S, -1.) );

  return grad ;
};

/*** R
library(nloptr)
library(purrr)
library(mvtnorm)
n <- 100; p <- 5; d <- 1
corr <- toeplitz(0.5^(1:p - 1))
sigma2 <- 2; rho <- 0.25
Sigma <- sigma2 * (rho * corr + (1 - rho) * diag(1, p, p))
Z <- rmvnorm(n, sigma = corr)
O <- matrix(log(round(runif(n, 10, 100))),n,p)
Y <- matrix(rpois(n* p, exp(O + Z)), n, p)
X <- matrix(1,n,d)
W <- matrix(rbinom(n * p, 1 , prob = 1/(1 + exp(- X %*% runif(d, -4, -2)))), n, p)
Y[W == 1] <- 0

LOGREGs <- lapply(1:p, function(j) glm.fit(X, 1*(Y[,j] != 0), offset =  O[,j], family = binomial(link = "logit")) )
Pi <- LOGREGs %>% map(fitted) %>% do.call(cbind, .)
Theta0 <- do.call(rbind, lapply(LOGREGs, coefficients))

LMs   <- lapply(1:p, function(j) lm.fit(X, log(1 + Y[,j]), offset =  O[,j]) )
Theta <- t(do.call(rbind, lapply(LMs, coefficients)))
M     <- do.call(cbind, lapply(LMs, residuals))
Omega <- solve(cov(M))

## Theta0

param <- rep(0, p*d)
print(objective_Theta0(param, X, Pi))
print(grad_Theta0(param, X, Pi))

objective_Theta0_R<- function(.x) {
  objective_Theta0(.x, X = X, Pi = Pi)
}

grad_Theta0_R<- function(.x) {
  grad_Theta0(.x, X = X, Pi = Pi)
}
check <- nloptr::check.derivatives(.x = param, objective_Theta0_R, grad_Theta0_R, check_derivatives_print = 'errors')


## M

param <- rep(0, n*p)
S <- sqrt(matrix(.1, n, p))

print(objective_M(param, Y, X, O, Pi, S, Theta, Omega))
print(grad_M(param, Y, X, O, Pi, S, Theta, Omega))

objective_M_R<- function(.x) {
  objective_M(.x, Y, X, O, Pi, S, Theta, Omega)
}

grad_M_R<- function(.x) {
  grad_M(.x, Y, X, O, Pi, S, Theta, Omega)
}
check <- nloptr::check.derivatives(.x = param, objective_M_R, grad_M_R, check_derivatives_print = 'errors')

## S
param <- rep(sqrt(0.1), n*p)

print(objective_S(param, O, M, Pi, Theta, Omega))
print(grad_S(param, O, M, Pi, Theta, Omega))

objective_S_R<- function(.x) {
  objective_S(.x, O, M, Pi, Theta, Omega)
}

grad_S_R<- function(.x) {
  grad_S(.x, O, M, Pi, Theta, Omega)
}
check <- nloptr::check.derivatives(.x = param, objective_S_R, grad_S_R, check_derivatives_print = 'errors')

*/
