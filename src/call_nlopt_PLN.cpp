#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "utils_optim.h"

double fn_optim_PLN(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++; // increase number of iterations

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::mat Omega = n * inv_sympd(M.t()*M  + diagmat(sum(S, 0)));

  double objective = accu(A - dat->Y % Z - .5*log(S)) - .5*n*real(log_det(Omega)) ; // + dat->KY ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise(M * Omega + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(Omega).t() + A - pow(S, -1)));

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}

double fn_optim_PLN_weighted(const std::vector<double> &x, std::vector<double> &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++; // increase number of iterations

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  double w_bar = accu(dat->w) ;

  arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % dat->w) + diagmat(sum(S.each_col() % dat->w, 0)));
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(diagmat(dat->w) *(A - dat->Y % Z - .5*log(S)) ) - .5 * w_bar*real(log_det(Omega)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M * Omega + A - dat->Y)) ;
  arma::vec grd_S     = vectorise(.5 * (dat->w * diagvec(Omega).t() + diagmat(dat->w) * A - diagmat(dat->w) * pow(S,-1) ) );

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}

arma::vec lower_bound_PLN_weighted(const std::vector<double> &x, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  double w_bar = accu(dat->w) ;

  arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % dat->w) + diagmat(sum(S.each_col() % dat->w, 0)));
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

 arma::vec J_i = sum(dat->Y % Z - A + .5*log(S), 1) + .5 * real(log_det(Omega)) - dat->KYi -
  .5 * (arma::diagvec(M * Omega * M.t()) + S * diagvec(Omega)) ;

  return J_i;
}

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
  nlopt::opt opt = initNLOPT(par.n_elem, options) ;

  // Perform the optimization
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);
  if (Rcpp::as<bool>(options["weighted"])) {
    opt.set_min_objective(fn_optim_PLN_weighted, &my_optim_data);
  } else {
    opt.set_min_objective(fn_optim_PLN, &my_optim_data);
  }

  nlopt::result status = opt.optimize(x_optimized, f_optimized);

  arma::vec lower_bound_terms = lower_bound_PLN_weighted(x_optimized, &my_optim_data) ;

  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("objective" ) = f_optimized + my_optim_data.KY ,
      Rcpp::Named("solution"  ) = x_optimized,
      Rcpp::Named("iterations") = my_optim_data.iterations,
      Rcpp::Named("loglik_obs") = lower_bound_terms
    );
}
