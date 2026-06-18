#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "builtin_optim_zipln.h"

// ---------------------------------------------------------------------------------------
// Builtin Newton VE-step for ZIPLN — all covariance structures (Omega and B fixed,
// joint (M, ψ, R) optimized). Each exported function delegates to the
// builtin_optimize_vestep_zipln_wrapper<Traits> template (builtin_optim_zipln.h), which
// extracts everything from (data, params, config).

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_zipln_full(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return builtin_optimize_vestep_zipln_wrapper<FullCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_zipln_diagonal(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return builtin_optimize_vestep_zipln_wrapper<DiagonalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_zipln_spherical(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return builtin_optimize_vestep_zipln_wrapper<SphericalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_zipln_fixed(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return builtin_optimize_vestep_zipln_wrapper<FixedCovTraits>(data, params, config);
}
