#include "gradients.h"

double fn_optim_PLN_full(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat S2 = S % S ;
  arma::mat A = exp (Z + .5 * S2) ;

  arma::mat Omega = dat->w_bar * inv_sympd(M.t() * (M.each_col() % dat->w) + diagmat((dat->w).t() * S2));

  double objective = accu((dat->w).t()*(A - dat->Y % Z - .5 * log(S2))) - .5*dat->w_bar*real(log_det(Omega)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M * Omega + A - dat->Y)) ;
  arma::vec grd_S     = vectorise(diagmat(dat->w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S,-1)) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S) ) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::vec S(&x[p*(d+n)], n);
  arma::vec S2 = S % S;
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S2) ;

  double sigma2 = arma::as_scalar(accu(M % (M.each_col() % dat->w)) / (dat->w_bar * p) + accu(dat->w % S2)/dat->w_bar);

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z))  - p*accu(.5 * dat->w % log(S2)) + .5 * dat->w_bar*p*log(sigma2) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M/sigma2 + A - dat->Y)) ;
  arma::vec grd_S     = dat->w % ( S % sum(A,1) -  p * pow(S, -1) - p * S/sigma2);

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}


double fn_optim_PLN_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  // recurrent variables
   int n = dat->n, p = dat->p, d = dat->d ;
  double w_bar = accu(dat->w) ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  arma::mat S2 = S % S;
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S2) ;

  arma::rowvec diag_Sigma = sum( M % (M.each_col() % dat->w) + (S2.each_col() % dat->w), 0) / w_bar;

  // objective function
  double objective = accu(diagmat(dat->w) *(A - dat->Y % Z - .5 * log(S2)) ) + .5 * w_bar* accu(log(diag_Sigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ( (M.each_row() / diag_Sigma) + A - dat->Y)) ;
  arma::vec grd_S     = vectorise(diagmat(dat->w) * (S.each_row() % pow(diag_Sigma, -1) + S % A - pow(S,-1)) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

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
  arma::mat S2 = S % S ;

  arma::mat Z = dat->O + dat->X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S2 * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_B     = vectorise((diagmat(dat->w) * (A - dat->Y)).t() * M + (A.t() * (S2.each_col() % dat->w)) % B) ;
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ((A-dat->Y) * B + M)) ;
  arma::vec grd_S     = vectorise(diagmat(dat->w) * (S - 1/S + A * (B%B) % S ));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S))) ;

  double objective = accu(diagmat(dat->w) * (A - dat->Y % Z)) + .5 * accu(diagmat(dat->w) * (M % M + S2 - log(S2) - 1)) ;

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
  arma::mat S2 = S % S ;

  arma::mat nSigma = M.t() * (M.each_col() % dat->w) + diagmat((dat->w).t() * S2);
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S2) ;

  double objective = accu((dat->w).t() * (A - dat->Y % Z - .5*log(S2))) - trace(dat->Omega*nSigma) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * (M * dat->Omega + A - dat->Y));
  arma::vec grd_S     = vectorise(diagmat(dat->w) * (S.each_row() % diagvec(dat->Omega).t() + S % A - pow(S,-1)) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M),grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_VEstep_PLN_full(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat S2 = S % S ;
  arma::mat Z = dat->O + dat->X * dat->Theta.t() + M;
  arma::mat A = exp (Z + .5 * S2) ;

  arma::mat mu = dat->X * dat->Theta.t() ;
  arma::mat nSigma =  M.t() * diagmat(dat->w) * M + diagmat(sum(S2.each_col() % dat->w, 0)) ;

  double objective = accu((dat->w).t()*(A - dat->Y % Z - .5*log(S2))) + .5*trace(dat->Omega*nSigma);

  arma::vec grd_M = vectorise(diagmat(dat->w) * ( M * dat->Omega + A - dat->Y)) ;
  arma::vec grd_S = vectorise(diagmat(dat->w) * (S.each_row() % diagvec(dat->Omega).t() + S % A - pow(S,-1)) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_VEstep_PLN_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  arma::mat M(&x[0]  , n,p);
  arma::vec S(&x[n*p], n);
  arma::vec S2 = S % S;

  double n_sigma2 = arma::as_scalar(dot(dat->w, sum(pow(M, 2), 1) + p * S)) ;
  double omega2 = arma::as_scalar(dat->Omega(0,0)) ;

  arma::mat Z = dat->O + dat->X * dat->Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S2) ;

  double objective = accu((dat->w).t() * (A - dat->Y % Z)) -.5 * p*dot(dat->w, log(S2)) +.5*n_sigma2*omega2 ;

  arma::vec grd_M = vectorise(diagmat(dat->w) * ( M * omega2 + A - dat->Y)) ;
  arma::vec grd_S = dat->w % ( S % sum(A,1) -  p * pow(S, -1) - p * S*omega2);

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;

  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_VEstep_PLN_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p ;

  // recurrent variables
  arma::mat M(&x[0]  , n,p);
  arma::mat S(&x[n*p], n,p);
  arma::mat Z = dat->O + dat->X * dat->Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;
  arma::mat S2 = S % S;

  arma::vec omega2 = arma::diagvec(dat->Omega);

  // objective function
  double objective = accu((dat->w).t()*(A - dat->Y % Z - .5*log(S2))) + .5 * as_scalar((dat->w).t() * (pow(M, 2) + S2) * omega2)  ;

  // gradients
  arma::vec grd_M = vectorise(diagmat(dat->w) * ( M * dat->Omega + A - dat->Y) ) ;
  arma::vec grd_S = vectorise(diagmat(dat->w) * (S.each_row() % omega2 + S % A - pow(S,-1)) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(grd_M, grd_S)) ;
  for (unsigned int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}
