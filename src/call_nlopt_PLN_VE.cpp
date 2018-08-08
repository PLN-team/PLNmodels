#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils_optim.h"

double fn_optim_VEstep_PLN(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++; // increase number of iterations

  int n = dat->Y.n_rows, p = dat->Y.n_cols ;

  arma::mat M(&x[0]    , n,p);
  arma::mat S(&x[n*p]  , n,p);

  arma::mat Z = dat->O + dat->X * dat->Theta.t() + M;
  arma::mat A = exp (Z + .5 * S);
  // 0.5 tr(\Omega M'M) + 0.5 tr(\bar{S} \Omega)
  double prior = .5*accu(dat->Omega % (M.t() * M)) + .5*dot(arma::ones(n).t() * S, diagvec(dat->Omega)) ;
  // J(M, S, \Theta, \Omega, Y, X, O)
  double objective = accu(A - dat->Y % Z - .5*log(S)) + prior - .5*n* dat->log_det_Omega + dat->KY ;

  arma::vec grd_M     = vectorise(M * dat->Omega + A-dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(dat->Omega).t() + A - 1/S));

  if (!grad.empty()) {
    grad = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  }

  return objective;
}

// [[Rcpp::export]]
Rcpp::List optimization_VEstep_PLN(
    arma::vec par,
    const arma::mat Y,
    const arma::mat X,
    const arma::mat O,
    const arma::mat Theta,
    const arma::mat Sigma,
    Rcpp::List options) {

  // Create optim data structure
  const arma::mat Omega = inv_sympd(Sigma);
  const double log_det_Omega = real(log_det(Omega));
  optim_data my_optim_data(Y, X, O, Theta, Omega, log_det_Omega);

  // Initialize the NLOPT optimizer
  nlopt::opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);
  opt.set_min_objective(fn_optim_VEstep_PLN, &my_optim_data);
  nlopt::result status = opt.optimize(x_optimized, f_optimized);

  // Format the output
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat M(&x_optimized[0]    , n,p);
  arma::mat S(&x_optimized[n*p]  , n,p);

  // Compute element-wise log-likelihood
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S);
  // sum(., 1) = rowSums(.)
  arma::vec loglik = arma::sum(-A + Y % Z + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)) - logfact(Y), 1) +
    .5*log_det_Omega*arma::ones(n);

  return Rcpp::List::create(
    Rcpp::Named("status"    ) = (int) status,
    Rcpp::Named("objective" ) = f_optimized,
    Rcpp::Named("M"         ) = M,
    Rcpp::Named("S"         ) = S,
    Rcpp::Named("loglik"    ) = loglik,
    Rcpp::Named("iterations") = my_optim_data.iterations
  );
}

