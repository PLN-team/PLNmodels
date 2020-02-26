#include "gradients.h"

double fn_optim_PLN(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z + .5 * S) ;

  arma::mat mu = dat->X * dat->Theta ;
  arma::mat nSigma =  (M-mu).t() * diagmat(dat->w) * (M-mu) + diagmat(sum(S.each_col() % dat->w, 0)) ;

  double objective = accu((dat->w).t()*(A - dat->Y % Z - .5*log(S))) + .5*trace(dat->Omega*nSigma) ;

  arma::vec grd_M = vectorise(diagmat(dat->w) * ( (M - mu) * dat->Omega + A - dat->Y)) ;
  arma::vec grd_S = vectorise(.5 * (dat->w * diagvec(dat->Omega).t() + diagmat(dat->w) * (A - pow(S,-1)) ) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::vec S(&x[n*p], n);

  arma::mat mu = dat->X * dat->Theta ;
  double n_sigma2 = arma::as_scalar(dot(dat->w, sum(pow(M - mu, 2), 1) + p * S)) ;
  double omega2 = arma::as_scalar(dat->Omega(0,0)) ;

  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;

  double objective = accu((dat->w).t() * (A - dat->Y % Z)) -.5 * p*dot(dat->w, log(S)) +.5*n_sigma2*omega2 ;

  arma::vec grd_M = vectorise(diagmat(dat->w) * ( (M - mu) * omega2 + A - dat->Y)) ;
  arma::vec grd_S = .5 * dat->w % (sum(A,1) +  p * omega2 - p * pow(S, -1));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}


double fn_optim_PLN_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  // recurrent variables
  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z + .5 * S) ;

  arma::mat mu = dat->X * dat->Theta ;
  arma::vec omega2 = arma::diagvec(dat->Omega);

  // objective function
  double objective = accu((dat->w).t()*(A - dat->Y % Z - .5*log(S))) + .5 * as_scalar((dat->w).t() * (pow(M - mu, 2) + S) * omega2)  ;

  // gradients
  arma::vec grd_M = vectorise(diagmat(dat->w) * ( (M - mu) * dat->Omega + A - dat->Y) ) ;
  arma::vec grd_S = vectorise(.5 * (dat->w * omega2.t() + diagmat(dat->w) * (A - pow(S,-1)) ) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_rank(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d, q = dat->q ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat B(&x[p*d]        , p,q) ;
  arma::mat M(&x[p*(d+q)]    , n,q) ;
  arma::mat S(&x[p*(d+q)+n*q], n,q) ;

  arma::mat Z = dat->O + dat->X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise((A-dat->Y).t() * dat->X);
  arma::vec grd_B     = vectorise((A-dat->Y).t() * M + (A.t() * S) % B) ;
  arma::vec grd_M     = vectorise((A-dat->Y) * B + M) ;
  arma::vec grd_S     = .5 * vectorise(1 - 1/S + A * (B%B) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S))) ;

  double objective = accu(A - dat->Y % Z) + .5 * accu(M % M + S - log(S) - 1) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_rank(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d, q = dat->q ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat B(&x[p*d]        , p,q) ;
  arma::mat M(&x[p*(d+q)]    , n,q) ;
  arma::mat S(&x[p*(d+q)+n*q], n,q) ;

  arma::mat Z = dat->O + dat->X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_B     = vectorise((diagmat(dat->w) * (A - dat->Y)).t() * M + (A.t() * (S.each_col() % dat->w)) % B) ;
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ((A-dat->Y) * B + M)) ;
  arma::vec grd_S     = .5 * vectorise(diagmat(dat->w) * (1 - 1/S + A * (B%B) ));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S))) ;

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z)) + .5 * accu(diagmat(dat->w) * (M % M + S - log(S) - 1)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_sparse(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat     M(&x[p*d]    , n,p) ;
  arma::mat     S(&x[p*(d+n)], n,p) ;

  arma::mat nSigma = M.t() * M ; nSigma.diag() += sum(S, 0);
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - dat->Y % Z - .5*log(S)) -.5*(n*dat->log_det_Omega + n*p - trace(dat->Omega*nSigma)) ;

  arma::vec grd_Theta = vectorise((A - dat->Y).t() * dat->X);
  arma::vec grd_M     = vectorise(M * dat->Omega + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(dat->Omega).t() + A - 1/S));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_sparse(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]      , p,d) ;
  arma::mat     M(&x[p*d]    , n,p) ;
  arma::mat     S(&x[p*(d+n)], n,p) ;
  double w_bar = accu(dat->w) ;

  arma::mat nSigma = M.t() * (M.each_col() % dat->w) + diagmat(sum(S.each_col() % dat->w, 0));
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z - .5*log(S))) -.5*(w_bar*dat->log_det_Omega + w_bar*p - trace(dat->Omega*nSigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M * dat->Omega + A - dat->Y));
  arma::vec grd_S     = vectorise(.5 * (dat->w *  diagvec(dat->Omega).t() + diagmat(dat->w) * A - diagmat(dat->w) * pow(S,-1)));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

