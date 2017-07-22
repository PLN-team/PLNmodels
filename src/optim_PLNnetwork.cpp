// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
Rcpp::List fn_optim_PLNnetwork_Cpp(const arma::vec par,
                               double log_detOmega,
                               const arma::mat Omega,
                               const arma::mat Y,
                               const arma::mat X,
                               const arma::mat O,
                               double KY) {

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta = par.subvec(0      , p*d      -1) ; Theta.reshape(p,d) ;
  arma::mat M     = par.subvec(p*d    , p*(n+d) - 1) ; M.reshape(n,p) ;
  arma::mat S     = par.subvec(p*(n+d), p*(2*n+d)-1) ; S.reshape(n,p) ;

  arma::mat S_bar(sum(S, 0)) ;
  arma::mat MtM = M.t() * M ;

  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double logP_Z = .5 * (n * log_detOmega - dot(diagvec(Omega), S_bar) - trace(Omega * MtM))  ;
  double objective = accu(A - Y % Z - .5 * log(S) - .5) - logP_Z + KY ;

  arma::vec grd_Theta = vectorise((A-Y).t() * X);
  arma::vec grd_M     = vectorise(M * Omega + A-Y) ;
  arma::vec grd_S     = vectorise(.5 * (ones(n) * diagvec(Omega).t() + A - 1/S));

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_M),grd_S) ;


  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}
