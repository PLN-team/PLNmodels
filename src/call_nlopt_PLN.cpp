#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "utils_optim.h"
#include "gradients_PLN.h"

// -----------------------------------------------
// FULLY PARAMETRIZED COVARIANCE MODEL

// [[Rcpp::export]]
Rcpp::List optimization_PLN(
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, w);

  // Initialize the NLOPT optimizer
  nlopt::opt optimizer = initNLOPT(par.n_elem, options) ;
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  // Perform the optimization
  if (Rcpp::as<bool>(options["weighted"])) {
    optimizer.set_min_objective(fn_optim_PLN_weighted, &my_optim_data);
  } else {
    optimizer.set_min_objective(fn_optim_PLN, &my_optim_data);
  }
  nlopt::result status = optimizer.optimize(x_optimized, f_optimized);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols;
  arma::mat Theta(&x_optimized[0]  , p,d);
  arma::mat M(&x_optimized[p*d]    , n,p);
  arma::mat S(&x_optimized[p*(d+n)], n,p);
  arma::mat Sigma = (M.t() * (M.each_col() % w) + diagmat(sum(S.each_col() % w, 0))) / accu(w) ;

  // Compute element-wise log-likelihood (works both for weighted/unweigthed loglikelihood)
  arma::mat Omega = inv_sympd(Sigma);
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::vec loglik = sum(Y % Z - A + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)), 1) + .5 * real(log_det(Omega)) - logfact(Y) + .5 * p;

  // Output returned to R
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized + my_optim_data.KY ,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("Theta" )     = Theta,
      Rcpp::Named("Sigma" )     = Sigma,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("iterations") = my_optim_data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}

// -----------------------------------------------
// DIAGONAL COVARIANCE MODEL

// [[Rcpp::export]]
Rcpp::List optimization_PLN_diagonal(
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, w);

  // Initialize the NLOPT optimizer
  nlopt::opt optimizer = initNLOPT(par.n_elem, options) ;
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  // Perform the optimization
  if (Rcpp::as<bool>(options["weighted"])) {
    optimizer.set_min_objective(fn_optim_PLN_weighted_diagonal, &my_optim_data);
  } else {
    optimizer.set_min_objective(fn_optim_PLN_diagonal, &my_optim_data);
  }
  nlopt::result status = optimizer.optimize(x_optimized, f_optimized);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols;
  arma::mat Theta(&x_optimized[0]  , p,d);
  arma::mat M(&x_optimized[p*d]    , n,p);
  arma::mat S(&x_optimized[p*(d+n)], n,p);
  arma::vec diag_Sigma = sum( M % (M.each_col() % w) + (S.each_col() % w), 0) / accu(w);
  arma::mat Sigma = diagmat(diag_Sigma) ;

  // Compute element-wise log-likelihood (works both for weighted/unweigthed loglikelihood)
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::vec loglik = sum(Y % Z - A + .5*log(S) - .5*( (M.each_row() / diag_Sigma) % M + (S.each_row() / diag_Sigma) ), 1) - .5 * sum(log(diag_Sigma)) - logfact(Y) + .5 * p;

  // Output returned to R
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized + my_optim_data.KY ,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("Theta" )     = Theta,
      Rcpp::Named("Sigma" )     = Sigma,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("iterations") = my_optim_data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}

// -----------------------------------------------
// SPHERICAL COVARIANCE

// [[Rcpp::export]]
Rcpp::List optimization_PLN_spherical(
    arma::vec par,
    const arma::mat & Y,
    const arma::mat & X,
    const arma::mat & O,
    const arma::vec & w,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, w);

  // Initialize the NLOPT optimizer
  nlopt::opt optimizer = initNLOPT(par.n_elem, options) ;
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  // Perform the optimization
  if (Rcpp::as<bool>(options["weighted"])) {
    optimizer.set_min_objective(fn_optim_PLN_weighted_spherical, &my_optim_data);
  } else {
    optimizer.set_min_objective(fn_optim_PLN_spherical, &my_optim_data);
  }
  nlopt::result status = optimizer.optimize(x_optimized, f_optimized);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols;
  arma::mat Theta(&x_optimized[0]  , p,d);
  arma::mat M(&x_optimized[p*d]    , n,p);
  arma::mat S(&x_optimized[p*(d+n)], n,1);
  double sigma2 = arma::as_scalar(accu(M % M) / (n * p) + accu(S)/n);
  arma::mat Sigma = arma::eye(p,p) * sigma2;

  // Compute element-wise log-likelihood (without weigths)
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;
  arma::vec loglik = sum(Y % Z - A, 1) + .5*p*log(S) - (.5 / sigma2) * (sum(M % M, 1) + p * S) - .5 * p * log(sigma2) - logfact(Y) + .5 * p;

  // Output returned to R
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized + my_optim_data.KY ,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("Theta" )     = Theta,
      Rcpp::Named("Sigma" )     = Sigma,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("iterations") = my_optim_data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}
