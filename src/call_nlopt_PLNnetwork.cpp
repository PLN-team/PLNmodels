#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "gradients_PLN.h"

// [[Rcpp::export]]
Rcpp::List optimization_PLNnetwork (
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::mat & Omega,
    const double log_det_Omega,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, Omega, log_det_Omega);

  // Initialize the NLOPT optimizer
  nlopt::opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);
  opt.set_min_objective(fn_optim_PLN_sparse, &my_optim_data);
  nlopt::result status = opt.optimize(x_optimized, f_optimized);

  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized  + my_optim_data.KY,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("iterations") = my_optim_data.iterations
    );
}
