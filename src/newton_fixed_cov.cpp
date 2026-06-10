#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl_alt.h"  // newton_optimize_alt_impl — the homemade E-step

// ---------------------------------------------------------------------------------------
// Fixed inverse covariance (Omega provided externally) PLN — homemade Newton optimizer

// [[Rcpp::export]]
Rcpp::List newton_optimize_fixed(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S, Omega)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S     = Rcpp::as<arma::mat>(params["S"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);

    const int    maxiter = config.containsElementNamed("maxeval")     ? Rcpp::as<int>(config["maxeval"])        : 200;
    const double ftol    = config.containsElementNamed("ftol_rel")    ? Rcpp::as<double>(config["ftol_rel"])    : 1e-8;
    const int    max_em  = config.containsElementNamed("maxit_em") ? Rcpp::as<int>(config["maxit_em"])    : 50;
    const double em_tol  = config.containsElementNamed("ftol_em")     ? Rcpp::as<double>(config["ftol_em"])     : 1e-8;

    FixedCovTraits::State state(Omega);
    return newton_optimize_alt_impl<FixedCovTraits>(Y, X, O, w, B, M, S, state, maxiter, ftol, max_em, em_tol);
}
