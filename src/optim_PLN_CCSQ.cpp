#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include <nlopt.hpp>

using namespace Rcpp;

typedef std::vector<double> stdvec;

struct optim_data {
    arma::mat Y   ;
    arma::mat X   ;
    arma::mat O   ;
    double KY     ;
    int iterations;
};

double fn_optim_PLN(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  // increase number of iterations
  dat->iterations++;

  int n = dat->Y.n_rows, p = dat->Y.n_cols, d = dat->X.n_cols ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Omega = n * inv_sympd(M.t()*M + diagmat(sum(S, 0)));
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - dat->Y % Z - .5*log(S)) - .5*n*real(log_det(Omega)) + dat->KY ;

  arma::vec grd_Theta = vectorise((A-dat->Y).t() * dat->X);
  arma::vec grd_M     = vectorise(M * Omega + A-dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(Omega).t() + A - 1/S));

  if (!grad.empty()) {
    grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;
  }

  return objective;
}

//' @export
// [[Rcpp::export]]
Rcpp::List optim_PLN_MMA(arma::vec par,
                         const arma::mat Y,
                         const arma::mat X,
                         const arma::mat O,
                         double         KY,
                         Rcpp::List control) {

  // Problem dimensions
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;
  int n_param = (2 * n + d) * p;

  // Create data structure
  optim_data my_optim_data;
  my_optim_data.Y  = Y  ;
  my_optim_data.X  = X  ;
  my_optim_data.O  = O  ;
  my_optim_data.KY = KY ; // compute this internally
  my_optim_data.iterations = 0 ; // compute this internally

  // prepare nlopt options
  double ftol_rel = as<double>(control["ftol_rel"]);
  double ftol_abs = as<double>(control["ftol_abs"]);
  double xtol_rel = as<double>(control["xtol_rel"]);
  double xtol_abs = as<double>(control["xtol_abs"]);
  int    maxeval  = as<int>   (control["maxeval"]);
  double lbvar    = as<double>(control["lbvar"]);
  stdvec lower_bound = arma::conv_to<stdvec>::from(arma::join_cols(R_NegInf * arma::ones(p*(d+n)), lbvar * arma::ones(n * p)));
  const stdvec xtol_abs_v = arma::conv_to<stdvec>::from(arma::join_cols(arma::zeros(p*(d+n)), xtol_abs * arma::ones(n * p)));

  arma::colvec test1 = arma::conv_to<arma::colvec>::from(lower_bound) ;
  arma::colvec test2 = arma::conv_to<arma::colvec>::from(xtol_abs_v) ;

  // Set nlopt options
  nlopt::opt opt(nlopt::LD_CCSAQ, n_param);

  opt.set_lower_bounds(lower_bound);
  opt.set_xtol_abs(xtol_abs_v);
  opt.set_xtol_rel(xtol_rel);
  opt.set_ftol_abs(ftol_abs);
  opt.set_ftol_rel(ftol_rel);
  opt.set_maxeval(maxeval);

  // Initialize the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  // Perform the optimization
  opt.set_min_objective(fn_optim_PLN, &my_optim_data);
  nlopt::result res = opt.optimize(x_optimized, f_optimized);
  int status  = (int) res ;

  return List::create(Named("status"    ) = status      ,
                      Named("objective" ) = f_optimized ,
                      Named("solution"  ) = x_optimized,
                      Named("iterations") = my_optim_data.iterations);
}
