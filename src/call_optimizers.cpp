#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

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
