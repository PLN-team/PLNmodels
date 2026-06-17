#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "builtin_optim_zipln.h"

// ---------------------------------------------------------------------------------------
// Builtin Newton VE-step for ZIPLN — all covariance structures (Omega and B fixed,
// joint (M, ψ, R) optimized). Each exported function extracts params/config (mirroring
// builtin_vestep_pln_impl's wrappers) and delegates to ve_step_zipln_newton_impl.

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_full(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    FullCovTraits::State state(Omega);
    return ve_step_zipln_newton_impl<FullCovTraits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_diagonal(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    DiagonalCovTraits::State state(Omega);
    return ve_step_zipln_newton_impl<DiagonalCovTraits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_spherical(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    SphericalCovTraits::State state(Omega);
    return ve_step_zipln_newton_impl<SphericalCovTraits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_fixed(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    FixedCovTraits::State state(Omega);
    return ve_step_zipln_newton_impl<FixedCovTraits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol);
}
