#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_pln.h"

// ─────────────────────────────────────────────────────────────────────────────
// ZIPLN VE-step math, reusing the CovTraitsBase<Traits> machinery from
// covariance_pln.h (Full/Diagonal/Spherical/Fixed) instead of duplicating a
// second hand-rolled set of formulas per covariance structure.
//
// Substitution that makes this possible: with A_eff = (1-R) ⊙ A (R fixed at
// the call site, e.g. for one Newton step or one inner nlopt solve) and
// w = ones(n) (ZIPLN carries no per-row weights, unlike PLN), the ZIPLN
// VE-step objective/gradient/Newton-step formulas are *identical* to
// CovTraitsBase::objective / ::compute_joint_step_MS with A replaced by
// A_eff — verified term-by-term. The only extra ingredient (the (1-R) split
// of the Poisson data term) is captured by zipln_vestep_obj_grad below for
// the nlopt path; the Newton path calls CovTraitsBase directly.
//
// The apparent mismatch where the ZIPLN ELBO weights *both* Y⊙Z and A by
// (1-R) while CovTraitsBase only reweights A is not a discrepancy: R is
// zero wherever Y > 0 by construction (R = σ(...) ⊙ 1{Y=0}), so Y⊙Z and
// (1-R)⊙Y⊙Z coincide identically.
// ─────────────────────────────────────────────────────────────────────────────

// Quantities derived from (Y, Pi) that stay constant through a VE-step call —
// computed once and shared by both backends (builtin_optim_zipln.h, nlopt_optim_zipln.h).
struct ZiplnRContext {
    arma::mat logit_Pi;
    arma::mat Y_zero;  // 1 where Y == 0, else 0 — R is restricted to these entries by construction
    ZiplnRContext(const arma::mat & Pi, const arma::mat & Y)
        : logit_Pi(logit(Pi)), Y_zero(arma::conv_to<arma::mat>::from(Y < 0.5)) {}
};

// Exact conditional optimum of R given (A, Pi): σ(A + logit(Pi)) where Y = 0, else 0.
inline arma::mat zipln_update_R(const arma::mat & A, const ZiplnRContext & ctx) {
    return (1.0 / (1.0 + arma::exp(-(A + ctx.logit_Pi)))) % ctx.Y_zero;
}

// VE-step objective + gradient for fixed R (A_eff = (1-R) ⊙ A passed in by the
// caller). Mirrors CovTraitsBase's private vestep_core, but takes A_eff as an
// external parameter instead of recomputing A = exp(Z + 0.5 S²) internally,
// since the (1-R) factor is not part of the State/Traits abstraction.
template <typename Traits>
inline double zipln_vestep_obj_grad(
    const arma::mat & M_res, const arma::mat & Z,
    const arma::mat & S2,   const arma::mat & logS2,
    const arma::mat & A_eff, const typename Traits::State & s,
    const arma::mat & Y,    const arma::vec & w,
    arma::mat & gM, arma::mat & gPS)
{
    const arma::mat MO = Traits::times_Omega(M_res, s);
    gM  = MO + A_eff - Y;                                        gM.each_col()  %= w;
    gPS = 0.5 * (Traits::diag_scale(S2, s) + S2 % A_eff - 1.0);  gPS.each_col() %= w;
    return arma::accu(w.t() * (A_eff - Y % Z - 0.5 * logS2))
         + CovTraitsBase<Traits>::penalty_M(MO, M_res, w) + Traits::penalty_S(S2, s, w);
}
