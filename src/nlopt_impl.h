#pragma once
#include <RcppArmadillo.h>

// Shared objective/gradient implementations for the three NLOPT covariance variants.
// Defined inline here so nlopt_full_cov.cpp, nlopt_diag_cov.cpp, nlopt_spherical.cpp and
// nlopt_fixed_cov.cpp can all include this header instead of each defining a private static copy.
//
// All three functions optimize over (M_full, ψ=log(S²)) and return grad_psi = d/dψ.

// Full covariance: Z = O + M_full supplied by caller.
// Returns the ELBO objective (negated); fills grad_M and grad_psi.
inline double full_cov_obj_grad_impl(
    const arma::mat & M_res, const arma::mat & Z, const arma::mat & logS2,
    const arma::mat & Omega, const arma::vec & Omega_diag,
    const arma::mat & Y,    const arma::vec & w,
    arma::mat & grad_M, arma::mat & grad_psi
) {
    const arma::mat S2 = arma::exp(logS2);
    const arma::mat A  = arma::exp(Z + 0.5 * S2);
    const arma::mat MO = M_res * Omega;
    grad_M   = MO + A - Y;                                              grad_M.each_col()   %= w;
    grad_psi = 0.5 * (S2.each_row() % Omega_diag.t() + S2 % A - 1.0); grad_psi.each_col() %= w;
    return accu(w.t() * (A - Y % Z - 0.5 * logS2))
         + 0.5 * (accu(MO % (M_res.each_col() % w)) + dot(Omega_diag, (w.t() * S2).t()));
}

// Diagonal covariance: inv_sigma2 = row vector of precisions (profiled or fixed);
// penalty = KL covariance term pre-computed by the caller (differs between E-step and vestep).
inline double diag_cov_obj_grad_impl(
    const arma::mat & M_res, const arma::mat & Z,
    const arma::mat & S2,   const arma::mat & logS2,
    const arma::rowvec & inv_sigma2, double penalty,
    const arma::mat & Y, const arma::vec & w,
    arma::mat & grad_M, arma::mat & grad_psi
) {
    const arma::mat A = arma::exp(Z + 0.5 * S2);
    grad_M   = M_res.each_row() % inv_sigma2 + A - Y;             grad_M.each_col()   %= w;
    grad_psi = 0.5 * (S2.each_row() % inv_sigma2 + S2 % A - 1.0); grad_psi.each_col() %= w;
    return accu(w.t() * (A - Y % Z - 0.5 * logS2)) + penalty;
}

// Spherical covariance: inv_sigma2 = scalar precision; penalty pre-computed by caller.
inline double spherical_cov_obj_grad_impl(
    const arma::mat & M_res, const arma::mat & Z,
    const arma::mat & S2,   const arma::mat & logS2,
    double inv_sigma2, double penalty,
    const arma::mat & Y, const arma::vec & w,
    arma::mat & grad_M, arma::mat & grad_psi
) {
    const arma::mat A = arma::exp(Z + 0.5 * S2);
    grad_M   = M_res * inv_sigma2 + A - Y;              grad_M.each_col()   %= w;
    grad_psi = 0.5 * (S2 * inv_sigma2 + S2 % A - 1.0); grad_psi.each_col() %= w;
    return accu(w.t() * (A - Y % Z - 0.5 * logS2)) + penalty;
}
