#include "optimizers.h"

// ---------------------------------------------------------------------------
// ABSTRACT CLASS OPTIMIZER_PLN
//
// COMMON TO PLN WITH SPHERICAL, DIAGONAL AND FULLY PARAMETRIZED COVARIANCE

// CONSTRUCTOR
optimizer_PLN::optimizer_PLN(
          arma::vec par,
          const arma::mat & Y,
          const arma::mat & X,
          const arma::mat & O,
          const arma::vec & w,
          Rcpp::List options) {

  // overload the data structure
  data = optim_data(Y, X, O, w) ;

  // problem dimension
  n = Y.n_rows ;
  p = Y.n_cols ;
  d = X.n_cols ;

  // Initialize NLOPT
  fn_optim   = NULL ;
  fn_VEstep  = NULL ;
  optimizer  = initNLOPT(par.n_elem, options)   ;
  parameter  = arma::conv_to<stdvec>::from(par) ;
}

// METHOD TO OPTIMIZE THE MAIN CRITERION
void optimizer_PLN::optimize()  {
  double objective ; // value of objective function at optimum

  nlopt_set_min_objective(optimizer, fn_optim, &data);
  status = nlopt_optimize(optimizer, &parameter[0], &objective) ;
  nlopt_destroy(optimizer);
}

// METHOD TO PERFORM A VE-STEP IN PLN
void optimizer_PLN::VEstep(const arma::mat & Theta, const arma::mat & Omega)  {
  double objective ; // value of objective function at optimum

  data.Theta = Theta ;
  data.Omega = Omega ;

  nlopt_set_min_objective(optimizer, fn_VEstep, &data);
  status = nlopt_optimize(optimizer, &parameter[0], &objective) ;
  nlopt_destroy(optimizer);
}

// CREATE THE RCPP::LIST FOR R
Rcpp::List optimizer_PLN::get_output() {
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("Theta"     ) = Theta,
      Rcpp::Named("Sigma"     ) = Sigma,
      Rcpp::Named("Omega"     ) = Omega,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("A"         ) = A,
      Rcpp::Named("Z"         ) = Z,
      Rcpp::Named("iterations") = data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}

// CREATE THE RCPP::LIST FOR R
Rcpp::List optimizer_PLN::get_var_par() {
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("iterations") = data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}

// ---------------------------------------------------------------------------
// CHILD CLASS WITH SPHERICAL COVARIANCE

optimizer_PLN_spherical::optimizer_PLN_spherical(
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {
  fn_optim  = &fn_optim_PLN_spherical  ;
  fn_VEstep = &fn_VEstep_PLN_spherical ;
}

void optimizer_PLN_spherical::export_output() {

  // variational parameters
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,1);
  arma::vec S2 = S % S ;

  // regression parameters
  Theta = arma::mat(&parameter[0]  , p,d);

  // variance parameters
  double sigma2 = arma::as_scalar(dot(data.w, sum(pow(M, 2), 1) + p * S2)) / (p * data.w_bar) ;
  Sigma = arma::eye(p,p) * sigma2 ;
  Omega = arma::eye(p,p) * pow(sigma2, -1) ;

  // element-wise log-likelihood
  Z = data.O + data.X * Theta.t() + M;
  A = exp(Z.each_col() + .5 * S2) ;
  loglik = sum(data.Y % Z - A - .5* pow(M, 2) / sigma2, 1) - .5 * p*S2/sigma2 + .5 *p*log(S2/sigma2) + data.Ki ;

}

void optimizer_PLN_spherical::export_var_par() {

  // variational parameters
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,1);
  arma::vec S2 = S % S ;

  double omega2 = arma::as_scalar(data.Omega(0,0)) ;

  // element-wise log-likelihood
  Z = data.O + data.X * Theta.t() + M;
  A = exp(Z.each_col() + .5 * S2) ;
  loglik = sum(data.Y % Z - A - .5* pow(M, 2) * omega2, 1) - .5*p*S2*omega2 + .5 *p*log(S2*omega2) + data.Ki ;
}

// ---------------------------------------------------------------------------
// CHILD CLASS WITH DIAGONAL COVARIANCE

optimizer_PLN_diagonal::optimizer_PLN_diagonal (
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {
  fn_optim  = &fn_optim_PLN_diagonal ;
  fn_VEstep = &fn_VEstep_PLN_diagonal ;
}

void optimizer_PLN_diagonal::export_output() {

  // variational parameters
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  arma::mat S2 = S % S ;

  // regression parameters
  Theta = arma::mat(&parameter[0]  , p,d);

  // variance parameters
  arma::rowvec sigma2 = (data.w).t() * (pow(M, 2) + S2) / data.w_bar;
  arma::vec omega2 = pow(sigma2.t(), -1) ;
  Sigma = diagmat(sigma2) ;
  Omega = diagmat(omega2) ;

  //element-wise log-likelihood
  Z = data.O + data.X * Theta.t() + M;
  A = exp (Z + .5 * S2) ;
  loglik = sum(data.Y % Z - A + .5 * log(S2), 1) - .5 * (pow(M, 2) + S2) * omega2 + .5 * sum(log(omega2)) + data.Ki ;
}

void optimizer_PLN_diagonal::export_var_par () {

  // variational parameters
  M = arma::mat(&parameter[0]  , n,p);
  S = arma::mat(&parameter[n*p], n,p);
  arma::vec omega2 = data.Omega.diag() ;
  arma::mat S2 = S % S ;

  //element-wise log-likelihood
  Z = data.O + data.X * Theta.t() + M;
  A = exp (Z + .5 * S2) ;
  loglik = sum(data.Y % Z - A + .5*log(S2), 1) - .5 * (pow(M, 2) + S2) * omega2 + .5 * sum(log(omega2)) + data.Ki ;
}


// ---------------------------------------------------------------------------
// CHILD CLASS WITH FULLY PARAMETRIZED COVARIANCE

optimizer_PLN_full::optimizer_PLN_full (
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {
  fn_optim  = &fn_optim_PLN_full ;
  fn_VEstep = &fn_VEstep_PLN_full ;
}

void optimizer_PLN_full::export_output () {

  // variational parameters
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  arma::mat S2 = S % S ;

  // regression parameters
  Theta = arma::mat(&parameter[0]  , p,d);

  // variance parameters
  Sigma = (M.t() * (M.each_col() % data.w) + diagmat(sum(S2.each_col() % data.w, 0))) / accu(data.w) ;
  Omega = inv_sympd(Sigma);

  // element-wise log-likelihood
  Z = data.O + data.X * Theta.t() + M ;
  A = exp (Z + .5 * S2) ;
  loglik = sum(data.Y % Z - A + .5* log(S2) - .5*( (M * Omega) % M + S2 * diagmat(Omega)), 1) + .5 * real(log_det(Omega)) + data.Ki ;
}

void optimizer_PLN_full::export_var_par () {

  // variational parameters
  M = arma::mat(&parameter[0]  , n,p);
  S = arma::mat(&parameter[n*p], n,p);
  arma::mat S2 = S % S;

  // element-wise log-likelihood
  Z = data.O + data.X * data.Theta.t() + M    ;
  A = exp (Z + .5 * S2) ;
  loglik = sum(data.Y % Z - A + .5*log(S2) - .5*( (M * data.Omega) % M + S * diagmat(data.Omega)), 1) + .5 * real(log_det(data.Omega)) + data.Ki ;
}

// ---------------------------------------------------------------------------
// CHILD CLASS WITH RANK-CONSTRAINED COVARIANCE
optimizer_PLN_rank::optimizer_PLN_rank (
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  const int rank,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {

  fn_optim = &fn_optim_PLN_rank ;

  // initialize the rank
  q = rank ;

  // complete the data structure
  data = optim_data(Y, X, O, w) ;
  data.q = rank ;
}

void optimizer_PLN_rank::export_output () {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]          , p,d);
  B     = arma::mat(&parameter[p*d]        , p,q) ;
  M     = arma::mat(&parameter[p*(d+q)]    , n,q);
  S     = arma::mat(&parameter[p*(d+q)+n*q], n,q);
  Z     = data.O + data.X * Theta.t() + M * B.t();
  arma::mat S2 = S % S ;
  Sigma = B * (M.t() * (M.each_col() % data.w) + diagmat(sum(S2.each_col() % data.w, 0))) * B.t() / accu(data.w) ;

  // element-wise log-likelihood
  A = exp (Z + .5 * S2 * (B % B).t() ) ;
  loglik = arma::sum(data.Y % Z - A, 1) - .5 * sum(M % M + S2 - log(S2) - 1, 1) + data.Ki;
}

// override mother's method for getting output
Rcpp::List optimizer_PLN_rank::get_output() {
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("Theta"     ) = Theta,
      Rcpp::Named("Sigma"     ) = Sigma,
      Rcpp::Named("B"         ) = B    ,
      Rcpp::Named("A"         ) = A    ,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S    ,
      Rcpp::Named("Z"         ) = Z,
      Rcpp::Named("iterations") = data.iterations,
      Rcpp::Named("loglik"    ) = loglik
    );
}

// ---------------------------------------------------------------------------
// CHILD CLASS WITH SPARSE INVERSE COVARIANCE

optimizer_PLN_sparse::optimizer_PLN_sparse (
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  const arma::mat & Omega,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {

  fn_optim = &fn_optim_PLN_sparse ;

  // complete the data structure
  data.Omega = Omega ;
  data.log_det_Omega = (real(log_det(Omega))) ;
}

void optimizer_PLN_sparse::export_output () {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]  , p,d);
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  arma::mat S2 = S % S ;
  Z = data.O + data.X * Theta.t() + M;
  A = exp (Z + .5 * S2) ;

  Sigma = (M.t() * (M.each_col() % data.w) + diagmat(data.w.t() * S2) )/ data.w_bar ;

  // element-wise log-likelihood
  loglik = sum(data.Y % Z - A - .5*( (M * data.Omega) % M - log(S2) + S2 * diagmat(data.Omega)), 1) + .5 * data.log_det_Omega  + data.Ki ;
}
