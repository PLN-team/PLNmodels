#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

// [[Rcpp::export]]
double objective_blocks(arma::vec param, const arma::mat &Y, const arma::mat&X, const arma::mat&O, const arma::mat&Tau) {
  arma::uword p = Y.n_cols;
  arma::uword d = X.n_cols;
  arma::uword n = Y.n_rows;
  arma::uword q = Tau.n_rows;

  arma::mat B(&param[0], d, p, true);
  arma::mat M(&param[p*d], n, q, true);
  arma::mat S(&param[p*d + n*q], n, q, true);

  arma::vec w = arma::ones(n) ;
  double w_bar = accu(w) ;
  arma::mat S2 = S % S;
  arma::mat mu = O + X * B ;
  arma::mat A1 = trunc_exp(M + .5 * S2) ;
  arma::mat A2 = trunc_exp(mu) ;
  arma::mat A = A2 % (A1 * Tau) ;
  arma::mat A_tau = ((A2 * Tau.t()) % A1) ;
  arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
  double objective = accu(w.t() * (A - Y % (mu + M * Tau))) - 0.5 * accu(w.t() * log(S2)) - 0.5 * w_bar * real(log_det(Omega));

  return objective ;
};

// [[Rcpp::export]]
arma::vec grad_blocks(arma::vec param, const arma::mat &Y, const arma::mat&X, const arma::mat&O, const arma::mat&Tau) {
  arma::uword p = Y.n_cols;
  arma::uword d = X.n_cols;
  arma::uword n = Y.n_rows;
  arma::uword q = Tau.n_rows;

  arma::mat B(&param[0], d, p, true);
  arma::mat M(&param[p*d], n, q, true);
  arma::mat S(&param[p*d + n*q], n, q, true);

  arma::vec w = arma::ones(n) ;
  double w_bar = accu(w) ;
  arma::mat S2 = S % S;
  arma::mat mu = O + X * B ;
  arma::mat A1 = trunc_exp(M + .5 * S2) ;
  arma::mat A2 = trunc_exp(mu) ;
  arma::mat A = A2 % (A1 * Tau) ;
  arma::mat A_tau = ((A2 * Tau.t()) % A1) ;
  arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));

  arma::vec B_grad = arma::vectorise((X.each_col() % w).t() * (A - Y));
  arma::vec M_grad =  arma::vectorise(diagmat(w) * (M * Omega + A_tau - Y * Tau.t()));
  arma::vec S_grad =  arma::vectorise(diagmat(w) * (S.each_row() % diagvec(Omega).t() + A_tau % S - pow(S, -1)));

  arma::vec grad = join_cols(B_grad, M_grad, S_grad);

  return grad ;
};


/*** R
library(nloptr)
data("oaks")
Y <- oaks$Abundance
X <- cbind(rep(1, nrow(Y)))
O <- log(oaks$Offset)

data("trichoptera")
Y <- as.matrix(trichoptera$Abundance)
X <- cbind(rep(1, nrow(Y)))
O <- matrix(0, n, p)

n <- nrow(Y); p <- ncol(Y); d <- ncol(X) ; q <- 5 ;

## Initialisation
fits <- lm.fit(X, log((1 + Y)/exp(O)))
B <- matrix(fits$coefficients, d, p)
M <- matrix(fits$residuals, n, p)
cl0 <- cutree(hclust(as.dist(1 - cor(M))), q)
Tau <- PLNmodels:::as_indicator(cl0) %>% t()
M <- M %*% t(Tau)
param <- c(B, M, matrix(0.1, n, q))

print(objective_blocks(param, Y, X, O, Tau))
print(grad_blocks(param, Y, X, O, Tau))

objective_blocks_R<- function(.x) {
  objective_blocks(.x, Y = Y, X = X, O = O, Tau = Tau)
}

grad_blocks_R<- function(.x) {
  grad_blocks(.x, Y = Y, X = X, O = O, Tau = Tau)
}

check <- nloptr::check.derivatives(.x = param, objective_blocks_R, grad_blocks_R)

plot(check$relative_error)
abline(v = c(p*d, p*d + n*q, p*d + 2*n*q))

*/