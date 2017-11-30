// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List fn_optim_PLN_fixpoint(arma::mat M,
                                 arma::mat S,
                                 arma::mat Theta,
                            const arma::mat Y,
                            const arma::mat X,
                            const arma::mat O,
                            double KY) {


  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat ProjX = X * inv_sympd( X.t() * X) * X.t();
  arma::mat ProjOrthX = diagmat(ones(n)) - X * inv_sympd( X.t() * X) * X.t();

  arma::mat A = exp (M + .5 * S), Omega, Sigma, Mtilde ;
  arma::vec objective(200);

  // update model parameters...
  Mtilde = ProjOrthX * (M - O);
  Sigma  = (Mtilde.t() * Mtilde + diagmat(sum(S, 0))) / n;
  Omega  = inv_sympd(Sigma);

  int i = 1;
  double dist = 2, eps = 0.01;

  objective[1] = accu(A - Y % M - .5*log(S)) - .5*n*real(log_det(Omega)) + KY ;
  std::cout << "\nitération: " << i << ", objective :" << objective[i] << endl;
  Sigma.print();

  while (dist > eps | i < 200) {
    i++ ;

    // update variational parameters...
    M = (Y-A) * Sigma + O + ProjX * (M-O);
    S = 1 / (A + ones(n) * diagvec(Omega).t()) ;

    // update model parameters...
    Mtilde = ProjOrthX * (M - O);
    Sigma  = (Mtilde.t() * Mtilde + diagmat(sum(S, 0))) / n;
    Sigma.print();

    Omega  = inv_sympd(Sigma);
    A = exp (M + .5 * S) ;

    objective[i] = accu(A - Y % M - .5*log(S)) - .5*n*real(log_det(Omega)) + KY ;
    std::cout << "\nitération: " << i << ", objective :" << objective[i] << "\n"<< endl;
//    Sigma.print();

//    A.print();
    dist = std::abs(objective(i)-objective(i-1))/objective(i);

  }
  Theta = (M-O).t() * X * inv_sympd(X.t() * X) ;
  return Rcpp::List::create(Rcpp::Named("M") = M,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("Theta") = Theta,
                            Rcpp::Named("Sigma") = Sigma,
                            Rcpp::Named("objective") = objective
                            );
}

/*** R
## Basic data simulation
n <- 10
p <- 5
d <- 1
rho <- 0.8
Sigma <- toeplitz(rho^(0:(p-1)))
Z <- matrix(rnorm(n*p), n,p) %*% chol(Sigma)
mu <- outer(rep(1,n), runif(p, 1, 2))
Y <- matrix(rpois(n*p, mu + exp(Z)), n,p)
X <- matrix(rep(1,n), n, 1)
O <- matrix(0,n,p)

## Set the generic
KY   <- sum(.logfactorial(Y))

PLN_out <- PLN(Y,X,O)
PLN_out
out <- fn_optim_PLN_fixpoint(PLN_out$var_par$M,
                             PLN_out$var_par$S,
                             PLN_out$model_par$Theta, Y, X, O, KY)
*/
