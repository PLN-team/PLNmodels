#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "builtin_optim_zipln.h"

// ---------------------------------------------------------------------------------------
// Builtin Newton VE-step for ZIPLN — all covariance structures (Omega and B fixed,
// joint (M, ψ, R) optimized). Each exported function delegates to the
// ve_step_zipln_newton_wrapper<Traits> template (builtin_optim_zipln.h), which extracts
// everything from (data, params, config).

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_full(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_newton_wrapper<FullCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_diagonal(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_newton_wrapper<DiagonalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_spherical(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_newton_wrapper<SphericalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton_fixed(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_newton_wrapper<FixedCovTraits>(data, params, config);
}
