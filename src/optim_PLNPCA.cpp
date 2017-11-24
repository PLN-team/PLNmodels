// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
Rcpp::List fn_optim_PLNPCA_Cpp(
    arma::vec par,
    int q,
    const arma::mat Y,
    const arma::mat X,
    const arma::mat O,
    double KY) {

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta(&par[0]          , p,d, false) ;
  arma::mat     B(&par[p*d]        , p,q, false) ;
  arma::mat     M(&par[p*(d+q)]    , n,q, false) ;
  arma::mat     S(&par[p*(d+q)+n*q], n,q, false) ;

  arma::mat Z = O + X * Theta.t() + M * B.t();
  arma::mat A = exp (Z + .5 * S * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise((A-Y).t() * X);
  arma::vec grd_B     = vectorise((A-Y).t() * M + (A.t() * S) % B) ;
  arma::vec grd_M     = vectorise((A-Y) * B + M) ;
  arma::vec grd_S     = .5 * vectorise(1 - 1/S + A * (B%B) );

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S)) ;

  double objective = accu(A - Y % Z) + .5 * accu(M%M + S - log(S) - 1) + KY ;

  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}
