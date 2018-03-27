#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include <nlopt.hpp>

using namespace Rcpp;

typedef std::vector<double> stdvec;

struct optim_data {
    arma::mat Y ;
    arma::mat X ;
    arma::mat O ;
    double KY;
};

double fn_optim_PLN(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);

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

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;

  return objective;
}


// static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
//     return (*reinterpret_cast<MyFunction*>(data))(x, grad);
// }

//' @export
// [[Rcpp::export]]
Rcpp::List optim_PLN_MMA(arma::vec par,
                         const arma::mat Y,
                         const arma::mat X,
                         const arma::mat O,
                         Rcpp::List control) {

  // Create data structure
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;
  optim_data my_optim_data;
  my_optim_data.Y  = Y  ;
  my_optim_data.X  = X  ;
  my_optim_data.O  = O  ;
  my_optim_data.KY = 1  ; // compute this internally


  int n_param = (2 * n + d) * p;

  nlopt::opt opt(nlopt::LD_MMA, n_param);

  opt.set_min_objective(fn_optim_PLN, NULL);

  opt.set_lower_bounds(lower_bound);

  opt.set_xtol_rel(1e-4);

  stdvec x = arma::conv_to<stdvec>::from(par);

  double minf;
  nlopt::result result = opt.optimize(x, minf);
  return List::create(Named("f") = minf, Named("x") = x);
}

/*** R
timesTwo(42)
*/
