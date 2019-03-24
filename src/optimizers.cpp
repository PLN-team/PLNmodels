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
          Rcpp::List options) : data(optim_data(Y, X, O, w)) {

  // problem dimension
  n = Y.n_rows ;
  p = Y.n_cols ;
  d = X.n_cols ;

  // Initialize NLOPT
  fn_optim = NULL ;
  optimizer = initNLOPT(par.n_elem, options)   ;
  parameter = arma::conv_to<stdvec>::from(par) ;
}

// FUNCTION THAT CALL NLOPT
void optimizer_PLN::optimize()  {
  double objective ; // value of objective function at optimum

  nlopt_set_min_objective(optimizer, fn_optim, &data);
  status = nlopt_optimize(optimizer, &parameter[0], &objective) ;
  nlopt_destroy(optimizer);
}

// CREATE THE RCPP::LIST FOR R
Rcpp::List optimizer_PLN::get_output() {
  return Rcpp::List::create(
      Rcpp::Named("status"    ) = (int) status,
      Rcpp::Named("Theta" )     = Theta,
      Rcpp::Named("Sigma" )     = Sigma,
      Rcpp::Named("M"         ) = M,
      Rcpp::Named("S"         ) = S,
      Rcpp::Named("A"         ) = A,
      Rcpp::Named("Z"         ) = Z,
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
  if (Rcpp::as<bool>(options["weighted"])) {
    fn_optim = &fn_optim_PLN_weighted_spherical ;
  } else {
    fn_optim = &fn_optim_PLN_spherical ;
  }
}

void optimizer_PLN_spherical::export_output() {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]  , p,d);
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,1);
  Z = data.O + data.X * Theta.t() + M;
  double sigma2 = arma::as_scalar(accu(M % M) / (accu(data.w) * p) + accu(S)/ accu(data.w)) ;
  Sigma = arma::eye(p,p) * sigma2;

  // element-wise log-likelihood
  A = exp (Z.each_col() + .5 * S) ;
  loglik = sum(data.Y % Z - A, 1) + .5*p*log(S) - (.5 / sigma2) * (sum(M % M, 1) + p * S) - .5 * p * log(sigma2) - logfact(data.Y) + .5 * p;
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
  if (Rcpp::as<bool>(options["weighted"])) {
    fn_optim = &fn_optim_PLN_weighted_diagonal ;
  } else {
    fn_optim = &fn_optim_PLN_diagonal ;
  }
}

void optimizer_PLN_diagonal::export_output() {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]  , p,d);
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  Z = data.O + data.X * Theta.t() + M;
  arma::rowvec diag_Sigma = sum( M % (M.each_col() % data.w) + (S.each_col() % data.w), 0) / accu(data.w);
  Sigma = diagmat(diag_Sigma) ;

  //element-wise log-likelihood
  A = exp (Z + .5 * S) ;
  loglik = sum(data.Y % Z - A + .5*log(S) - .5*( (M.each_row() / diag_Sigma) % M + (S.each_row() / diag_Sigma) ), 1) - .5 * accu(log(diag_Sigma)) - logfact(data.Y) + .5 * p;
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
  if (Rcpp::as<bool>(options["weighted"])) {
    fn_optim = &fn_optim_PLN_weighted ;
  } else {
    fn_optim = &fn_optim_PLN ;
  }
}

void optimizer_PLN_full::export_output () {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]  , p,d);
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  Z = data.O + data.X * Theta.t() + M    ;
  Sigma = (M.t() * (M.each_col() % data.w) + diagmat(sum(S.each_col() % data.w, 0))) / accu(data.w) ;

  // element-wise log-likelihood
  arma::mat Omega = inv_sympd(Sigma);
  A = exp (Z + .5 * S) ;
  loglik = sum(data.Y % Z - A + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)), 1) + .5 * real(log_det(Omega)) - logfact(data.Y) + .5 * p;
}

// ---------------------------------------------------------------------------
// CHILD CLASS WITH RANK-CONSTRAINED COVARIANCE

optimizer_PLN_rank::optimizer_PLN_rank (
  arma::vec par,
  const arma::mat & Y,
  const arma::mat & X,
  const arma::mat & O,
  const arma::vec & w,
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {
  if (Rcpp::as<bool>(options["weighted"])) {
    fn_optim = &fn_optim_PLN_weighted_rank ;
  } else {
    fn_optim = &fn_optim_PLN_rank ;
  }

  // initialize the rank
  q = Rcpp::as<int>(options["rank"]);

  // overload the data structure
  data = optim_data(Y, X, O, w, q) ;

}

void optimizer_PLN_rank::export_output () {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]          , p,d);
  B     = arma::mat(&parameter[p*d]        , p,q) ;
  M     = arma::mat(&parameter[p*(d+q)]    , n,q);
  S     = arma::mat(&parameter[p*(d+q)+n*q], n,q);
  Z     = data.O + data.X * Theta.t() + M * B.t();
  Sigma = B * (M.t()* M + diagmat(sum(S, 0)) ) * B.t() / n ;

  // element-wise log-likelihood
  A = exp (Z + .5 * S * (B % B).t() ) ;
  loglik = arma::sum(data.Y % Z - A, 1) - .5 * sum(M % M + S - log(S) - 1, 1) - logfact(data.Y);
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
      Rcpp::Named("S"         ) = S,
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
  Rcpp::List options
) : optimizer_PLN(par, Y, X, O, w, options) {
  if (Rcpp::as<bool>(options["weighted"])) {
    fn_optim = &fn_optim_PLN_weighted_sparse ;
  } else {
    fn_optim = &fn_optim_PLN_sparse ;
  }

  const arma::mat & Omega = Rcpp::as<arma::mat>(options["Omega"]);

  // overload the data structure
  data = optim_data(Y, X, O, w, Omega) ;

}

void optimizer_PLN_sparse::export_output () {

  // model and variational parameters
  Theta = arma::mat(&parameter[0]  , p,d);
  M = arma::mat(&parameter[p*d]    , n,p);
  S = arma::mat(&parameter[p*(d+n)], n,p);
  Z = data.O + data.X * Theta.t() + M;
  Sigma = (M.t() * (M.each_col() % data.w) + diagmat(sum(S.each_col() % data.w, 0))) / accu(data.w) ;

  // element-wise log-likelihood
  A = exp (Z + .5 * S) ;
  loglik = sum(data.Y % Z - A + .5*log(S) - .5*( (M * data.Omega) % M + S * diagmat(data.Omega)), 1) + .5 * data.log_det_Omega - logfact(data.Y) + .5 * p;
}
