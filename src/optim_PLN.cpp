// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
Rcpp::List fn_optim_PLN_par2_Cpp(arma::vec par,
                            const arma::mat Y,
                            const arma::mat ProjOrthX,
                            const arma::mat O,
                            double KY) {
  // Gradients checked numerically

  int n = Y.n_rows, p = Y.n_cols ;

  arma::mat M(&par[0]  , n,p, false) ;
  arma::mat S(&par[n*p], n,p, false) ;

   // SVD of Mtilde * S-1/2
   // arma::mat U, V;
   // arma::vec d ;
   // svd_econ(U, d, V, ProjOrthX * (M - O) * diagmat(1/sqrt(s)));

   arma::mat A = exp (M + .5 * S) ;
   arma::mat Mtilde = ProjOrthX * (M - O);

   arma::mat Omega  = n * inv_sympd( Mtilde.t() * Mtilde + diagmat(sum(S, 0)));

   double objective = accu(A - Y % M - .5*log(S)) - .5*n*real(log_det(Omega)) + KY ;
   arma::vec grd_M     = vectorise( Mtilde * Omega + A  - Y) ;
   arma::vec grd_S     = vectorise(.5 * (ones(n) * diagvec(Omega).t() + A - 1/S));

  // double objective = accu(A - Y % M - .5*log(S)) + .5* n * sum(log( s/n % (1+square(d)))) + KY ;
  //
  // arma::mat VtS = V.t() * diagmat(1/sqrt(s)) ;
  // arma::vec grd_M     = vectorise( n * U * diagmat(d/(1 + square(d))) * VtS + A - Y) ;
  //
  //  vec diagOmega = diagvec( n * VtS.st() * diagmat(1/(1+square(d))) * VtS ) ;
  //  arma::vec grd_S = vectorise(.5 * (ones(n) * diagOmega.t() + A - 1/S));

  arma::vec grad = join_vert(grd_M,grd_S) ;

  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}

//' @export
// [[Rcpp::export]]
Rcpp::List fn_optim_PLN_par1_Cpp( arma::vec par,
                                  const arma::mat Y,
                                  const arma::mat X,
                                  const arma::mat O,
                                  double KY) {
  // Gradients checked numerically

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta(&par[0]  , p,d, false);
  arma::mat M(&par[p*d]    , n,p, false);
  arma::mat S(&par[p*(d+n)], n,p, false);

  arma::mat Omega = n * inv_sympd(M.t()*M + diagmat(sum(S, 0)));
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - Y % Z - .5*log(S)) - .5*n*real(log_det(Omega)) + KY ;

  arma::vec grd_Theta = vectorise((A-Y).t() * X);
  arma::vec grd_M     = vectorise(M * Omega + A-Y) ;
  arma::vec grd_S     = vectorise(.5 * (ones(n) * diagvec(Omega).t() + A - 1/S));

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_M),grd_S) ;

  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}
