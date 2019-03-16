#include "gradients.h"

double fn_optim_PLN(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
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

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
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

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::vec S(&x[p*(d+n)], n);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z.each_col() + .5 * S) ;
  double sigma2 = arma::as_scalar(accu(M % M) / (n * p) + accu(S)/n);

  double objective = accu(A - dat->Y % Z)  - .5*p*accu(log(S)) + .5*n*p*std::log(sigma2) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise(M/sigma2 + A - dat->Y) ;
  arma::vec grd_S     = .5 * (sum(A,1) -  p * pow(S, -1) - p/sigma2);

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_spherical(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
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

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}


double fn_optim_PLN_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  arma::rowvec diag_Sigma = sum(M % M  + S, 0) / n;

  double objective = accu(A - dat->Y % Z - .5*log(S)) + .5*n*accu(log(diag_Sigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * dat->X);
  arma::vec grd_M     = vectorise((M.each_row() / diag_Sigma) + A - dat->Y) ;
  arma::vec grd_S     = vectorise(.5 * (arma::ones(n) * pow(diag_Sigma, -1) + A - pow(S, -1)));

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}

double fn_optim_PLN_weighted_diagonal(unsigned N, const double *x, double *grad, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);
  double w_bar = accu(dat->w) ;

  arma::rowvec diag_Sigma = sum( M % (M.each_col() % dat->w) + (S.each_col() % dat->w), 0) / w_bar;
  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(diagmat(dat->w) *(A - dat->Y % Z - .5*log(S)) ) + .5 * w_bar* accu(log(diag_Sigma)) ;

  arma::vec grd_Theta = vectorise(trans(A - dat->Y) * (dat->X.each_col() % dat->w));
  arma::vec grd_M     = vectorise(diagmat(dat->w) * ( (M.each_row() / diag_Sigma) + A - dat->Y)) ;
  arma::vec grd_S     = vectorise(.5 * (dat->w * pow(diag_Sigma, -1) + diagmat(dat->w) * A - diagmat(dat->w) * pow(S,-1) ) );

  stdvec grad_std = arma::conv_to<stdvec>::from(join_vert(join_vert(grd_Theta, grd_M), grd_S)) ;

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

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

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

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

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

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

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

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

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

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

  for (int i=0;i<N;i++) grad[i] = grad_std[i];

  return objective;
}
