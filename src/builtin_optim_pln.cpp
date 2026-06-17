#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "builtin_optim_pln.h"

// ---------------------------------------------------------------------------------------
// Builtin Newton optimizer for PLN — all covariance structures.
// Each exported function is a thin wrapper that builds the appropriate Traits::State
// and delegates to the generic builtin_optimize_pln_impl / builtin_vestep_pln_impl templates.

// ===== FULL COVARIANCE =====

// [[Rcpp::export]]
Rcpp::List builtin_optimize_full(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);
    const NewtonConfig cfg(config);
    const double w_bar = arma::accu(d.w);
    FullCovTraits::State state(M - d.X * B, S2, d.w, w_bar);
    return builtin_optimize_pln_impl<FullCovTraits>(d, B, M, S2, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_full(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2    = Rcpp::as<arma::mat>(params["S2"]);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    FullCovTraits::State state(Omega);
    return builtin_vestep_pln_impl<FullCovTraits>(d, M, S2, B, state, cfg.maxiter, cfg.ftol);
}

// ===== DIAGONAL COVARIANCE =====

// [[Rcpp::export]]
Rcpp::List builtin_optimize_diagonal(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);
    const NewtonConfig cfg(config);
    const double w_bar = arma::accu(d.w);
    DiagonalCovTraits::State state(M - d.X * B, S2, d.w, w_bar);
    return builtin_optimize_pln_impl<DiagonalCovTraits>(d, B, M, S2, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_diagonal(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2    = Rcpp::as<arma::mat>(params["S2"]);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    DiagonalCovTraits::State state(Omega);
    return builtin_vestep_pln_impl<DiagonalCovTraits>(d, M, S2, B, state, cfg.maxiter, cfg.ftol);
}

// ===== SPHERICAL COVARIANCE =====

// [[Rcpp::export]]
Rcpp::List builtin_optimize_spherical(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);
    const NewtonConfig cfg(config);
    const double w_bar = arma::accu(d.w);
    SphericalCovTraits::State state(M - d.X * B, S2, d.w, w_bar);
    return builtin_optimize_pln_impl<SphericalCovTraits>(d, B, M, S2, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_spherical(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2    = Rcpp::as<arma::mat>(params["S2"]);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    SphericalCovTraits::State state(Omega);
    return builtin_vestep_pln_impl<SphericalCovTraits>(d, M, S2, B, state, cfg.maxiter, cfg.ftol);
}

// ===== FIXED COVARIANCE (Omega provided externally, no VE step exported) =====

// [[Rcpp::export]]
Rcpp::List builtin_optimize_fixed(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData d(data);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2    = Rcpp::as<arma::mat>(params["S2"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    FixedCovTraits::State state(Omega);
    return builtin_optimize_pln_impl<FixedCovTraits>(d, B, M, S2, state, cfg.maxiter, cfg.ftol, cfg.max_em, cfg.em_tol);
}
