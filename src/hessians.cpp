#include "hessians.h"

void fn_precond_PLN(unsigned N, const double *x, const double *v, double *Hv, void *data) {

  optim_data *dat = (optim_data *) data;
  dat->iterations++;

  int n = dat->n, p = dat->p, d = dat->d ;

  arma::mat X_M(&v[0]  , n,p);
  arma::mat X_S(&v[p*n], n,p);

  arma::mat Theta(&x[0]  , p,d);
  arma::mat M(&x[p*d]    , n,p);
  arma::mat S(&x[p*(d+n)], n,p);

  arma::mat Z = dat->O + dat->X * Theta.t() + M;
  arma::mat Omega = n * inv_sympd(M.t()*M  + diagmat(sum(S, 0)));
  arma::mat AMS = exp (Z + .5 * S) % (X_M + .5 * X_S);

  arma::mat Hv_M =      AMS + X_M * Omega  ;
  arma::mat Hv_S = .5* (AMS + X_S /(S % S)) ;

  stdvec Hv_std = arma::conv_to<stdvec>::from(join_vert(Hv_M, Hv_S)) ;

  for (unsigned int i=0;i<N;i++) Hv[i] = Hv_std[i];

}
