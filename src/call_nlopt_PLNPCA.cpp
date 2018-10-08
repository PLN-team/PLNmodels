#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "utils_optim.h"

double fn_optim_PLNPCA(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d, q = dat->q ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat B(&x[p*d]        , p,q) ;
  arma::mat M(&x[p*(d+q)]    , n,q) ;
  arma::mat S(&x[p*(d+q)+n*q], n,q) ;

  arma::mat Z = dat->O + dat->X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise((A-dat->Y).t() * dat->X);
  arma::vec grd_B     = vectorise((A-dat->Y).t() * M + (A.t() * S) % B) ;
  arma::vec grd_M     = vectorise((A-dat->Y) * B + M) ;
  arma::vec grd_S     = .5 * vectorise(1 - 1/S + A * (B%B) );

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S))) ;

  double objective = accu(A - dat->Y % Z) + .5 * accu(M % M + S - log(S) - 1) ;

  return objective;
}

// [[Rcpp::export]]
Rcpp::List optimization_PLNPCA (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const int rank,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, rank);

  // Initialize the NLOPT optimizer
  nlopt::opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);
  opt.set_min_objective(fn_optim_PLNPCA, &my_optim_data);
  nlopt::result status = opt.optimize(x_optimized, f_optimized);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols, q = rank;
  arma::mat Theta(&x_optimized[0]      , p,d) ;
  arma::mat B(&x_optimized[p*d]        , p,q) ;
  arma::mat M(&x_optimized[p*(d+q)]    , n,q) ;
  arma::mat S(&x_optimized[p*(d+q)+n*q], n,q) ;
  arma::mat Sigma = B * (M.t()* M + diagmat(sum(S, 0)) ) * B.t() / n ;

  // Compute element-wise lower bound of the log-likelihood
  arma::mat Z = O + X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S * (B % B).t() ) ;
  arma::vec loglik = arma::sum(Y % Z - A, 1) - .5 * sum(M % M + S - log(S) - 1, 1) - logfact(Y);

  return Rcpp::List::create(
      Rcpp::Named("status") = (int) status,
      Rcpp::Named("loglik")    = loglik,
      Rcpp::Named("Theta" )    = Theta,
      Rcpp::Named("B"     )    = B,
      Rcpp::Named("Sigma" )    = Sigma,
      Rcpp::Named("M"     )    = M,
      Rcpp::Named("S"     )    = S,
      Rcpp::Named("iterations") = my_optim_data.iterations
    );
}
