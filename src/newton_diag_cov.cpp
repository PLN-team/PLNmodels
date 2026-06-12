#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"  // newton_optimize_impl + newton_vestep_impl

// ---------------------------------------------------------------------------------------
// Diagonal covariance PLN — homemade Newton optimizer (profiled B via envelope theorem)

// [[Rcpp::export]]
Rcpp::List newton_optimize_diagonal(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);

    const NewtonConfig cfg(config);

    const double w_bar = arma::accu(w);
    const arma::mat M_res_init = M - X * B;
    DiagonalCovTraits::State state(M_res_init, S2, w, w_bar);

    return newton_optimize_impl<DiagonalCovTraits>(Y, X, O, w, B, M, S2, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol);
}

// ---------------------------------------------------------------------------------------
// VE diagonal — coordinate-Newton (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List newton_optimize_vestep_diagonal(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p) fixed
    const arma::mat & Omega,   // (p,p) diagonal, fixed
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);

    const NewtonConfig cfg(config);
    DiagonalCovTraits::State state(Omega);
    return newton_vestep_impl<DiagonalCovTraits>(Y, X, O, w, M, S2, B, state, cfg.maxiter, cfg.ftol);
}
