#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "utils_optim.h"

double fn_optim_PLN(const stdvec &x, stdvec &grad, void *data) {

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

double fn_optim_PLN_weighted(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

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

double fn_optim_PLN_spherical(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++; // increase number of iterations

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::vec S(&x[p*(d+n)], n);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;
  double sigma2 = arma::as_scalar(accu(M % M) / (n * p) + accu(S)/n);

  double objective = accu(A - dat->Y % Z)  - .5*p*accu(log(S)) - .5*n*p*std::log(sigma2) ; // + dat->KY ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise(M/sigma2 + A - dat->Y) ;
  arma::vec grd_S     = .5 * (sum(A,1) -  p * pow(S, -1) - p/sigma2);

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
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
  nlopt::opt optimizer = initNLOPT(par.n_elem, options) ;
  double f_optimized ;
  stdvec x_optimized = arma::conv_to<stdvec>::from(par);

  // Perform the optimization
  if (Rcpp::as<bool>(options["weighted"])) {
    optimizer.set_min_objective(fn_optim_PLN_weighted, &my_optim_data);
  } else if (Rcpp::as<std::string>(options["covariance"]) == "spherical") {
    optimizer.set_min_objective(fn_optim_PLN_spherical, &my_optim_data);
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
  arma::vec loglik = sum(Y % Z - A + .5*log(S) - .5*( (M * Omega) % M + S * diagmat(Omega)), 1) +
    + .5 * real(log_det(Omega)) - logfact(Y) + my_optim_data.w * .5 * p;

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
