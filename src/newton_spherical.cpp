#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"  // newton_optimize_impl + newton_vestep_impl

// ---------------------------------------------------------------------------------------
// Spherical covariance PLN — homemade Newton optimizer (profiled B via envelope theorem)

// [[Rcpp::export]]
Rcpp::List newton_optimize_spherical(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const NewtonConfig cfg(config);

    const double w_bar = arma::accu(w);
    arma::mat S2 = S % S;
    const arma::mat M_res_init = M - X * B;
    SphericalCovTraits::State state(M_res_init, S2, w, w_bar);

    return newton_optimize_impl<SphericalCovTraits>(Y, X, O, w, B, M, S, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol, cfg.block_newton_thresh);
}

// ---------------------------------------------------------------------------------------
// VE spherical — coordinate-Newton (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List newton_optimize_vestep_spherical(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p) fixed
    const arma::mat & Omega,   // (p,p) diagonal scalar, fixed
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const NewtonConfig cfg(config);
    SphericalCovTraits::State state(Omega);
    return newton_vestep_impl<SphericalCovTraits>(Y, X, O, w, M, S, B, state, cfg.maxiter, cfg.ftol, cfg.block_newton_thresh);
}
