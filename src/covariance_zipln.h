#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_pln.h"

// ─────────────────────────────────────────────────────────────────────────────
// ZIPLN VE-step math, on top of CovTraitsBase<Traits> (covariance_pln.h).
//
// At fixed R and with unit weights, the objective, gradient and Newton step of
// the ZIPLN VE step are those of CovTraitsBase with A replaced by
// A_eff = (1-R) ⊙ A. The ZIPLN ELBO also weights Y⊙Z by (1-R), which changes
// nothing: R is zero wherever Y > 0. zipln_vestep_obj_grad below is the nlopt
// path; the Newton path calls CovTraitsBase directly.
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

// VE-step objective + gradient for fixed R: as CovTraitsBase's vestep_core,
// with A_eff = (1-R) ⊙ A passed in by the caller
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
