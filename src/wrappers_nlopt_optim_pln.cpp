#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "nlopt_optim_pln.h"

// ---------------------------------------------------------------------------------------
// Full covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const PlnData d(data);
    const auto init_B  = Rcpp::as<arma::mat>(params["B"]);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);

    // Initial Omega profiled from (M_res_init, S2_init); the M-step then re-profiles
    // it after every inner nlopt run (see nlopt_optimize_em_impl).
    const arma::mat M_res_init = init_M - d.X * init_B;
    FullCovTraits::State state(M_res_init, init_S2, d.w, arma::accu(d.w));

    return nlopt_optimize_em_impl<FullCovTraits>(data, init_M, init_S2, state, config);
}

// ---------------------------------------------------------------------------------------
// Full covariance PLN — profiled variant (no EM loop): Omega profiled at every eval.
// Recommended for small p where the per-eval O(np² + p³) cost is acceptable.
// For large p, prefer nlopt_optimize_full (EM with fixed-Omega inner nlopt).

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_full_profiled(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    return nlopt_joint_profiled_impl<FullCovTraits>(data, params, config);
}

// ---------------------------------------------------------------------------------------
// VE full covariance — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_full(
    const Rcpp::List & data,   // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S2, B, Omega)
    const Rcpp::List & config
) {
    return nlopt_vestep_impl<FullCovTraits>(data, params, config);
}

// ---------------------------------------------------------------------------------------
// Fixed covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_fixed(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);

    FixedCovTraits::State state(Omega);  // Omega given externally, never updated (has_em = false)
    return nlopt_optimize_em_impl<FixedCovTraits>(data, init_M, init_S2, state, config);
}

// ---------------------------------------------------------------------------------------
// Diagonal covariance PLN — nlopt/CCSAQ optimizer: B and sigma² profiled at every eval

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_diagonal(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    return nlopt_joint_profiled_impl<DiagonalCovTraits>(data, params, config);
}

// ---------------------------------------------------------------------------------------
// VE diagonal — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_diagonal(
    const Rcpp::List & data,   // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S2, B, Omega)
    const Rcpp::List & config
) {
    return nlopt_vestep_impl<DiagonalCovTraits>(data, params, config);
}

// ---------------------------------------------------------------------------------------
// Spherical covariance PLN — nlopt/CCSAQ optimizer: B and sigma² profiled at every eval

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_spherical(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    return nlopt_joint_profiled_impl<SphericalCovTraits>(data, params, config);
}

// ---------------------------------------------------------------------------------------
// VE spherical — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_spherical(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S2, B, Omega)
    const Rcpp::List & config
) {
    return nlopt_vestep_impl<SphericalCovTraits>(data, params, config);
}
