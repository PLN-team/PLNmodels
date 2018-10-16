#include "gradients_PLN.h"

double fn_optim_PLN(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::mat Omega = n * inv_sympd(M.t()*M  + diagmat(sum(S, 0)));

  double objective = accu(A - dat->Y % Z - .5*log(S)) - .5*n*real(log_det(Omega)) ;

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

  double objective = accu(A - dat->Y % Z)  - .5*p*accu(log(S)) + .5*n*p*std::log(sigma2) ; // + dat->KY ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise(M/sigma2 + A - dat->Y) ;
  arma::vec grd_S     = .5 * (sum(A,1) -  p * pow(S, -1) - p/sigma2);

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}

double fn_optim_PLN_weighted_spherical(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::vec S(&x[p*(d+n)], n);
  double w_bar = accu(dat->w) ;

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;
  double sigma2 = arma::as_scalar(accu(M % (M.each_col() % dat->w)) / (w_bar * p) + accu(dat->w % S)/w_bar);

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z))  - .5*p*accu(dat->w % log(S)) + .5*w_bar*p*log(sigma2) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M/sigma2 + A - dat->Y)) ;
  arma::vec grd_S     = dat->w % (.5 * (sum(A,1) -  p * pow(S, -1) - p/sigma2));

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}


double fn_optim_PLN_diagonal(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  arma::vec diag_Sigma = sum(M % M  + S, 0) / n;

  double objective = accu(A - dat->Y % Z - .5*log(S)) + .5*n*accu(log(diag_Sigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise((M.each_row() / diag_Sigma) + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * pow(diag_Sigma, -1).t() + A - pow(S, -1)));

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}

double fn_optim_PLN_weighted_diagonal(const stdvec &x, stdvec &grad, void *data) {

  optim_data *dat = reinterpret_cast<optim_data*>(data);
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  double w_bar = accu(dat->w) ;

  arma::vec diag_Sigma = sum( M % (M.each_col() % dat->w) + (S.each_col() % dat->w), 0) / w_bar;
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(diagmat(dat->w) *(A - dat->Y % Z - .5*log(S)) ) + .5 * w_bar* accu(log(diag_Sigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ( (M.each_row() / diag_Sigma) + A - dat->Y)) ;
  arma::vec grd_S     = vectorise(.5 * (dat->w * pow(diag_Sigma, -1).t() + diagmat(dat->w) * A - diagmat(dat->w) * pow(S,-1) ) );

  grad = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  return objective;
}
