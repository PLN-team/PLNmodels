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

  arma::mat nSigma = (M - mu).t()*(M-mu)  + diagmat(sum(S, 0));

  double objective = accu(A - dat->Y % Z - .5*log(S)) - .5*(n*dat->log_det_Omega - trace(dat->Omega*nSigma)) ;

  arma::vec grd_M     = vectorise( (M - mu) * dat->Omega + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(dat->Omega).t() + A - pow(S, -1))) ;

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::mat mu = dat->X * dat->XtWX_inv * dat->X.t() * (M.each_col() % dat->w) ;

  arma::mat Omega = dat->w_bar * inv_sympd( (M-mu).t() * diagmat(dat->w) * (M-mu) + diagmat(sum(S.each_col() % dat->w, 0))) ;

  double objective = accu(diagmat(dat->w) *(A - dat->Y % Z - .5*log(S)) ) - .5 * dat->w_bar * real(log_det(Omega)) ;

  arma::vec grd_M = vectorise(diagmat(dat->w) * ( (M - mu) * Omega + A - dat->Y)) ;
  arma::vec grd_S = vectorise(.5 * (dat->w * diagvec(Omega).t() + diagmat(dat->w) * (A - pow(S,-1)) ) );

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

  // regression parameters
  arma::mat Mu = dat->X * dat->XtWX_inv * dat->X.t() * M ;

  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;
  double sigma2 = arma::as_scalar(accu(M % M) / (n * p) + accu(S)/n);

  double objective = accu(A - dat->Y % Z) - .5*p*accu(log(S/sigma2)) ;

  arma::vec grd_M = vectorise((M - Mu)/sigma2 + A - dat->Y ) ;
  arma::vec grd_S = .5 * (sum(A,1) +  p/sigma2 - p * pow(S, -1));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::vec S(&x[n*p], n);

  // regression parameters
  arma::mat mu = dat->X * dat->XtWX_inv * dat->X.t() * (M.each_col() % dat->w) ;

  // variance parameters
  double sigma2 = arma::as_scalar(dot(dat->w, sum(pow(M - mu, 2), 1) + p * S)) / (p * dat->w_bar) ;

  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z)) -.5 *p*dot(dat->w, log(S/sigma2)) ;

  arma::vec grd_M     = vectorise(diagmat(dat->w) * ( (M - mu)/sigma2 + A - dat->Y)) ;
  arma::vec grd_S     = .5 * dat->w % (sum(A,1) +  p/sigma2 - p * pow(S, -1));

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

  // regression parameters
  arma::mat mu = dat->X * dat->XtWX_inv * dat->X.t() * M ;

  // variance parameters
  arma::rowvec d = sum(pow(M - mu, 2) + S, 0) / n;

  // objective function
  double objective = accu(A - dat->Y % Z - .5*log(S)) + .5 * n * accu(log(d)) ;

  // gradients
  arma::vec grd_M = vectorise( (M - mu) * diagmat(pow(d, -1)) + A - dat->Y ) ;
  arma::vec grd_S = vectorise(.5 * (arma::ones(n) * pow(d, -1) + A - pow(S,-1) ) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  // recurrent variables
  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat Z = dat->O + M;
  arma::mat A = exp (Z + .5 * S) ;

  // regression parameters
  arma::mat mu = dat->X * dat->XtWX_inv * dat->X.t() * (M.each_col() % dat->w) ;

  // variance parameters
  arma::rowvec d = (dat->w).t() * (pow(M - mu, 2) + S) / dat->w_bar;

  // objective function
  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z - .5*log(S))) + .5 * dat->w_bar * accu(log(d)) ;

  // gradients
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ( (M - mu) * diagmat(pow(d, -1)) + A - dat->Y) ) ;
  arma::vec grd_S     = vectorise(.5 * (dat->w * pow(d, -1) + diagmat(dat->w) * (A - pow(S,-1)) ) );

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

double fn_optim_VEstep_PLN(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++; // increase number of iterations

  int n = dat->Y.n_rows, p = dat->Y.n_cols ;

  arma::mat M(&x[0]    , n,p);
  arma::mat S(&x[n*p]  , n,p);

  arma::mat Z = dat->O + dat->X * dat->Theta.t() + M;
  arma::mat A = exp (Z + .5 * S);
  // 0.5 tr(\Omega M'M) + 0.5 tr(\bar{S} \Omega)
  double prior = .5*accu(dat->Omega % (M.t() * M)) + .5*dot(arma::ones(n).t() * S, diagvec(dat->Omega)) ;
  // J(M, S, \Theta, \Omega, Y, X, O)
  double objective = accu(A - dat->Y % Z - .5*log(S)) + prior - .5*n* dat->log_det_Omega ;

  arma::vec grd_M     = vectorise(M * dat->Omega + A-dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * diagvec(dat->Omega).t() + A - 1/S));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}
