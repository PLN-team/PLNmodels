#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "optimizer_PLN.h"

// [[Rcpp::export]]
Rcpp::List optimization_PLN (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Initialize the optimizer
  optimizer_PLN *myPLN ;

  // SPHERICAL COVARIANCE
  if (Rcpp::as<std::string>(options["covariance"]) == "spherical")
    myPLN = new optimizer_PLN_spherical(par, Y, X, O, w, options) ;

  // DIAGONAL COVARIANCE
  if (Rcpp::as<std::string>(options["covariance"]) == "diagonal")
    myPLN = new optimizer_PLN_diagonal(par, Y, X, O, w, options) ;

  // FULLY PARAMETRIZED COVARIANCE
  if (Rcpp::as<std::string>(options["covariance"]) == "full")
    myPLN = new optimizer_PLN_full(par, Y, X, O, w, options) ;

  // Perform the optimization
  myPLN->optimize() ;

  // Format the output
  myPLN->export_output() ;

  // Output returned to R
  return(myPLN->get_output());

}


