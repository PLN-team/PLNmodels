#include "RcppArmadillo.h"
#include "nloptrAPI.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "data_struct.h"
#include "nlopt_utils.h"
#include "optimizers.h"

// ---------------------------------------------------------------------------------------
// SPHERICAL COVARIANCE
//
//
// [[Rcpp::export]]
Rcpp::List optim_spherical (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialization
  optimizer_PLN_spherical myPLN = optimizer_PLN_spherical(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN.optimize() ;

  // Format the output
  myPLN.export_output() ;

  // Output returned to R
  return(myPLN.get_output());
}

// ---------------------------------------------------------------------------------------
// DIAGONAL COVARIANCE
//
//
// [[Rcpp::export]]
Rcpp::List optim_diagonal (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialize
  optimizer_PLN_diagonal myPLN = optimizer_PLN_diagonal(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN.optimize() ;

  // Format the output
  myPLN.export_output() ;

  // Output returned to R
  return(myPLN.get_output());
}

// ---------------------------------------------------------------------------------------
// FULLY PARAMETRIZED COVARIANCE
//
//
// [[Rcpp::export]]
Rcpp::List optim_full (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialize
  optimizer_PLN_full myPLN = optimizer_PLN_full(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN.optimize() ;

  // Format the output
  myPLN.export_output() ;

  // Output returned to R
  return(myPLN.get_output());
}

// ---------------------------------------------------------------------------------------
// RANK-CONSTRAINED COVARIANCE (PCA)
//
//
// [[Rcpp::export]]
Rcpp::List optim_rank (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialize
  optimizer_PLN_rank myPLN = optimizer_PLN_rank(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN.optimize() ;

  // Format the output
  myPLN.export_output() ;

  // Output returned to R
  return(myPLN.get_output());
}

// ---------------------------------------------------------------------------------------
// SPARSE INVERSE COVARIANCE (aka 'NETWORK')
//
//
// [[Rcpp::export]]
Rcpp::List optim_sparse (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialize
  optimizer_PLN_sparse myPLN = optimizer_PLN_sparse(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN.optimize() ;

  // Format the output
  myPLN.export_output() ;

  // Output returned to R
  return(myPLN.get_output());
}

// function to perform a single VE Step

// [[Rcpp::export]]
Rcpp::List optimization_VEstep_PLN(
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::mat & Theta,
    const arma::mat & Sigma,
    Rcpp::List options) {

  // Create optim data structure
  const arma::mat Omega = inv_sympd(Sigma);
  const double log_det_Omega = real(log_det(Omega));
  optim_data my_optim_data(Y, X, O, Theta, Omega, log_det_Omega);

  // Initialize the NLOPT optimizer
  nlopt_opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  nlopt_set_min_objective(opt, fn_optim_VEstep_PLN, &my_optim_data);
  nlopt_result status = nlopt_optimize(opt, &x_optimized[0], &f_optimized);
  nlopt_destroy(opt);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat M(&x_optimized[0]    , n,p);
  arma::mat S(&x_optimized[n*p]  , n,p);

  // Compute element-wise log-likelihood
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S);
  // sum(., 1) = rowSums(.)
  arma::vec loglik = arma::sum(-A + Y % Z + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)), 1) +
    .5*log_det_Omega - logfact(Y) + .5 * p;

  return Rcpp::List::create(
    Rcpp::Named("status"    ) = (int) status,
    Rcpp::Named("objective" ) = f_optimized + my_optim_data.KY ,
    Rcpp::Named("M"         ) = M,
    Rcpp::Named("S"         ) = S,
    Rcpp::Named("loglik"    ) = loglik,
    Rcpp::Named("iterations") = my_optim_data.iterations
  );
}

