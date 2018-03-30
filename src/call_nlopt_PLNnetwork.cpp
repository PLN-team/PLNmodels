#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils_optim.h"

double fn_optim_PLNnetwork(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++; // increase number of iterations

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat     M(&x[p*d]    , n,p) ;
  arma::mat     S(&x[p*(d+n)], n,p) ;

  arma::mat nSigma = M.t() * M ; nSigma.diag() += sum(S, 0);
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - dat->Y % Z - .5*log(S)) -.5*(n*dat->log_det_Omega + n*p - trace(dat->Omega*nSigma)) + dat->KY ;

  arma::vec grd_Theta = vectorise((A - dat->Y).t() * dat->X);
  arma::vec grd_M     = vectorise(M * dat->Omega + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(dat->Omega).t() + A - 1/S));

  if (!grad.empty()) {
    grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;
  }

  return objective;
}

//' @export
// [[Rcpp::export]]
Rcpp::List optimization_PLNnetwork (
    arma::vec par,
    const arma::mat Y,
    const arma::mat X,
    const arma::mat O,
    const arma::mat Omega,
    const double log_det_Omega,
    Rcpp::List options) {

  // Create data structure
  optim_data my_optim_data(Y, X, O, Omega, log_det_Omega);

  // Initialize the NLOPT optimizer
  nlopt::opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);
  opt.set_min_objective(fn_optim_PLNnetwork, &my_optim_data);
  nlopt::result status = opt.optimize(x_optimized, f_optimized);

  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized ,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("iterations") = my_optim_data.iterations
    );
}
