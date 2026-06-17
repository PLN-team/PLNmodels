#pragma once
#include <RcppArmadillo.h>
#include "utils.h"

// ─────────────────────────────────────────────────────────────────────────────
// Rank-constrained covariance — shared math for the joint (nlopt/builtin) and
// VE-step (nlopt/builtin) rank-constrained PLN optimizers.
//
// Unlike full/diagonal/spherical/fixed PLN, B has no closed-form profiling here:
// the KL term only involves the q-dim score M (standard normal prior, independent
// of B), and B appears only in the non-quadratic Poisson data term via
// Z = O + X*B + M*C'. So B stays in the optimized parameter vector for both
// backends, and there is exactly one covariance structure (rank-constrained) —
// no CRTP/trait abstraction is warranted here, just the shared objective/gradient
// math that both backends would otherwise duplicate.
//
// Variational parameter: ψ = log(S²) (unconstrained) instead of S2 (bounded > 0).
// ─────────────────────────────────────────────────────────────────────────────

// Joint objective + gradient (B, C, M, ψ all free).
// Xw = X.each_col() % w, precomputed once by the caller and reused across evals.
inline double rank_obj_grad(
    const PlnData & d, const arma::mat & Xw,
    const arma::mat & B, const arma::mat & C, const arma::mat & M, const arma::mat & psi,
    arma::mat & gB, arma::mat & gC, arma::mat & gM, arma::mat & gPS)
{
    const arma::mat C2  = C % C;
    const arma::mat S2  = arma::exp(psi);
    const arma::mat Z   = d.O + d.X * B + M * C.t();
    const arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    const arma::mat AmY = A - d.Y;

    double objective = arma::accu(d.w.t() * (A - d.Y % Z))
                      + 0.5 * arma::accu(d.w.t() * (M % M + S2 - psi - 1.));

    arma::mat AmYw = AmY; AmYw.each_col() %= d.w;
    gB  = Xw.t() * AmY;
    gC  = AmYw.t() * M + (A.t() * (S2.each_col() % d.w)) % C;
    gM  = AmY * C + M;                           gM.each_col()  %= d.w;
    gPS = 0.5 * (S2 % (1. + A * C2) - 1.);        gPS.each_col() %= d.w;

    return objective;
}

// VE-step objective + gradient (B, C fixed; M, ψ free).
// XB = X*B and C2 = C%C precomputed once by the caller.
inline double rank_vestep_obj_grad(
    const PlnData & d, const arma::mat & XB, const arma::mat & C, const arma::mat & C2,
    const arma::mat & M, const arma::mat & psi,
    arma::mat & gM, arma::mat & gPS)
{
    const arma::mat S2  = arma::exp(psi);
    const arma::mat Z   = d.O + XB + M * C.t();
    const arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    const arma::mat AmY = A - d.Y;

    double objective = arma::accu(d.w.t() * (A - d.Y % Z))
                      + 0.5 * arma::accu(d.w.t() * (M % M + S2 - psi - 1.));

    gM  = AmY * C + M;                    gM.each_col()  %= d.w;
    gPS = 0.5 * (S2 % (1. + A * C2) - 1.); gPS.each_col() %= d.w;

    return objective;
}

// Per-observation log-likelihood, shared by the joint and VE-step variants.
inline arma::vec rank_final_loglik(
    const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
    const arma::mat & M, const arma::mat & S2, const arma::mat & psi)
{
    return arma::sum(Y % Z - A, 1) - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(Y);
}

// Sigma/Omega derived from the converged (M, C, S2) — joint variant only
// (the VE-step keeps B, C, and hence Sigma/Omega, fixed).
inline Rcpp::List rank_output_cov(
    const arma::mat & M, const arma::mat & C, const arma::mat & S2,
    const arma::vec & w, double w_bar)
{
    const arma::mat nSig  = M.t() * (M.each_col() % w) + arma::diagmat(arma::sum(S2.each_col() % w, 0));
    const arma::mat Sigma = C * nSig * C.t() / w_bar;
    const arma::mat Omega = C * arma::inv_sympd(nSig / w_bar) * C.t();
    return Rcpp::List::create(Rcpp::Named("Sigma", Sigma), Rcpp::Named("Omega", Omega));
}
