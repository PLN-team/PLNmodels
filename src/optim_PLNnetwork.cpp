// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
Rcpp::List fn_optim_PLNnetwork_Cpp(
    arma::vec par,
    double log_detOmega,
    const arma::mat Omega,
    const arma::mat Y,
    const arma::mat X,
    const arma::mat O,
    double KY) {

  // Gradient checked - ok up to 1e-3
  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta(&par[0]      , p,d, false) ;
  arma::mat     M(&par[p*d]    , n,p, false) ;
  arma::mat     S(&par[p*(d+n)], n,p, false) ;

  arma::mat nSigma = M.t() * M ; nSigma.diag() += sum(S, 0);
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - Y % Z - .5*log(S)) -.5*(n*log_detOmega + n*p - trace(Omega*nSigma)) + KY ;

  arma::vec grd_Theta = vectorise((A-Y).t() * X);
  arma::vec grd_M     = vectorise(M * Omega + A-Y) ;
  arma::vec grd_S     = vectorise(.5 * (ones(n) * diagvec(Omega).t() + A - 1/S));

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_M),grd_S) ;

  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}
