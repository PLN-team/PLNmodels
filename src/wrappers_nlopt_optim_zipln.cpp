#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "nlopt_optim_zipln.h"

// ---------------------------------------------------------------------------------------
// nlopt VE-step for ZIPLN — all covariance structures (Omega and B fixed, joint
// (M, ψ, R) optimized). Each exported function is a thin wrapper that delegates to
// ve_step_zipln_nlopt_impl, which extracts everything from (data, params, config).

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_nlopt_full(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_nlopt_impl<FullCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_nlopt_diagonal(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_nlopt_impl<DiagonalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_nlopt_spherical(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_nlopt_impl<SphericalCovTraits>(data, params, config);
}

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_nlopt_fixed(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    return ve_step_zipln_nlopt_impl<FixedCovTraits>(data, params, config);
}
