#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"  // newton_optimize_impl — the homemade E-step

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

    const NewtonConfig cfg(config);

    FixedCovTraits::State state(Omega);
    return newton_optimize_impl<FixedCovTraits>(Y, X, O, w, B, M, S, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol, cfg.block_newton_thresh);
}
