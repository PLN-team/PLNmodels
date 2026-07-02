#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_plnpca.h"

// ─────────────────────────────────────────────────────────────────────────────
// Profiled trust-region Newton for rank-constrained PLN (PLNPCA), fixed rank q.
//
// Idea (validated in R prototypes): the variational block (M, ψ) is a CONCAVE
// inner problem at fixed (B, C); profile it out with a per-observation Newton
// VE-step, giving g(B,C) = max_{M,ψ} ELBO and, by the envelope theorem, ∇g for
// free (= ∂ELBO/∂(B,C) at the VE optimum). Then optimise g(B,C) with a
// saddle-aware trust-region Newton, whose reduced (Schur-complement) Hessian
//     H_red = L_θθ − L_θφ L_φφ⁻¹ L_φθ ,   θ = (B,C),  φ = (M,ψ)
// is applied matrix-free: L_φφ is block-diagonal per observation (2q×2q,
// analytic below); the cross/θθ terms use directional finite differences of the
// analytic gradient `rank_obj_grad`. Block-coordinate ascent gets stuck at
// saddles of this non-convex landscape; the negative-curvature-aware TR does not.
// ─────────────────────────────────────────────────────────────────────────────

namespace trnewton {

// Objective (minimisation convention, matching rank_obj_grad) at fixed (B,C).
inline double ve_obj(const PlnData & d, const arma::mat & XB, const arma::mat & C,
                     const arma::mat & C2, const arma::mat & M, const arma::mat & psi) {
    const arma::mat S2 = arma::exp(psi);
    const arma::mat Z  = d.O + XB + M * C.t();
    const arma::mat A  = arma::exp(arma::clamp(Z + 0.5 * S2 * C2.t(), -arma::datum::inf, 30.0));
    return arma::accu(d.w.t() * (A - d.Y % Z)) + 0.5 * arma::accu(d.w.t() * (M % M + S2 - psi - 1.));
}

// VE-step: minimise over (M, ψ) at fixed (B, C). Per-observation Newton on M
// (q×q solve, concave block) + exact fixed-point for S² = exp(ψ). Modifies M, ψ.
inline double ve_solve(const PlnData & d, const arma::mat & XB, const arma::mat & C,
                       const arma::mat & C2, arma::mat & M, arma::mat & psi,
                       int maxit = 50, double tol = 1e-10) {
    const arma::uword n = d.Y.n_rows, q = C.n_cols;
    const arma::mat Iq = arma::eye(q, q);
    double f0 = ve_obj(d, XB, C, C2, M, psi);
    if (!std::isfinite(f0)) { M.zeros(); psi.fill(std::log(0.1)); f0 = ve_obj(d, XB, C, C2, M, psi); }
    for (int it = 0; it < maxit; ++it) {
        arma::mat S2 = arma::exp(psi);
        const arma::mat Z = d.O + XB + M * C.t();
        arma::mat A = arma::exp(arma::clamp(Z + 0.5 * S2 * C2.t(), -arma::datum::inf, 30.0));
        // exact fixed-point for S² given M (∂/∂ψ = 0)
        S2  = 1.0 / (1.0 + A * C2);
        psi = arma::log(S2);
        A   = arma::exp(arma::clamp(Z + 0.5 * S2 * C2.t(), -arma::datum::inf, 30.0));
        // Newton step on M, per observation (weights cancel in the step)
        const arma::mat gM = (A - d.Y) * C + M;      // gradient of the (min) objective
        arma::mat dM(n, q);
        for (arma::uword i = 0; i < n; ++i) {
            arma::mat DC = C; DC.each_col() %= A.row(i).t();
            arma::mat H  = C.t() * DC + Iq;
            dM.row(i) = arma::solve(H, gM.row(i).t(),
                                    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx).t();
        }
        double t = 1.0, f1 = ve_obj(d, XB, C, C2, M - t * dM, psi);
        while ((!std::isfinite(f1) || f1 > f0 + 1e-12) && t > 1e-6) { t *= 0.5; f1 = ve_obj(d, XB, C, C2, M - t * dM, psi); }
        M -= t * dM;
        if (std::isfinite(f1) && std::abs(f1 - f0) <= tol * std::abs(f1)) { f0 = f1; break; }
        f0 = f1;
    }
    return f0;
}

// Analytic per-observation inner Hessian L_φφ (2q×2q, φ = (M,ψ)), returned as a
// cube of its inverses (one slice per observation) for the Schur solve.
//   H_MM = Cᵀdiag(A_i)C + I
//   H_Mψ = ½ (Cᵀdiag(A_i)C²)·diag(S²)
//   H_ψψ = ½ diag(S²(1+G)) + ¼ diag(S²)(C²ᵀdiag(A_i)C²)diag(S²),  G = C²ᵀA_i
inline arma::cube inner_blocks_inv(const PlnData & d, const arma::mat & XB, const arma::mat & C,
                                   const arma::mat & C2, const arma::mat & M, const arma::mat & psi) {
    const arma::uword n = d.Y.n_rows, q = C.n_cols;
    const arma::mat S2 = arma::exp(psi);
    const arma::mat Z  = d.O + XB + M * C.t();
    const arma::mat A  = arma::exp(arma::clamp(Z + 0.5 * S2 * C2.t(), -arma::datum::inf, 30.0));
    const arma::mat Iq = arma::eye(q, q), I2 = arma::eye(2 * q, 2 * q);
    arma::cube Hinv(2 * q, 2 * q, n);
    for (arma::uword i = 0; i < n; ++i) {
        const arma::rowvec ai = A.row(i);
        const arma::vec    s2 = S2.row(i).t();
        arma::mat DC  = C;  DC.each_col()  %= ai.t();
        arma::mat DC2 = C2; DC2.each_col() %= ai.t();
        const arma::mat P  = C.t()  * DC;    // Cᵀdiag(a)C
        arma::mat Q        = C.t()  * DC2;   // Cᵀdiag(a)C²  (Q_{kl}=Σ a C_jk C_jl²)
        const arma::mat Rm = C2.t() * DC2;   // C²ᵀdiag(a)C²
        const arma::vec G  = C2.t() * ai.t();
        arma::mat Hmm = P + Iq;
        arma::mat Hmp = Q; Hmp.each_row() %= s2.t(); Hmp *= 0.5;              // ½ Q diag(S²)
        arma::mat Hpp = 0.5 * arma::diagmat(s2 % (1.0 + G))
                      + 0.25 * (arma::diagmat(s2) * Rm * arma::diagmat(s2));
        arma::mat Hi(2 * q, 2 * q);
        Hi.submat(0, 0, q - 1, q - 1)         = Hmm;
        Hi.submat(0, q, q - 1, 2 * q - 1)     = Hmp;
        Hi.submat(q, 0, 2 * q - 1, q - 1)     = Hmp.t();
        Hi.submat(q, q, 2 * q - 1, 2 * q - 1) = Hpp;
        Hi *= d.w(i);
        Hinv.slice(i) = arma::inv_sympd(Hi + 1e-8 * I2);
    }
    return Hinv;
}

// Hessian of the joint objective applied to a direction (dB,dC,dM,dψ), returning
// the induced change of every gradient block. The only nonlinearity is A=exp(η),
// so everything reduces to dA = A ⊙ dη + product rule. Consistent term-by-term
// with `rank_obj_grad`. A, S2 are the values at the current (B,C,M,ψ).
inline void hess_dir(const PlnData & d, const arma::mat & Xw,
                     const arma::mat & C, const arma::mat & C2, const arma::mat & M,
                     const arma::mat & A, const arma::mat & S2,
                     const arma::mat & dB, const arma::mat & dC,
                     const arma::mat & dM, const arma::mat & dpsi,
                     arma::mat & dgB, arma::mat & dgC, arma::mat & dgM, arma::mat & dgPS) {
    const arma::vec & w = d.w;
    const arma::mat dS2  = S2 % dpsi;           // n×q
    const arma::mat dC2  = 2.0 * (C % dC);      // p×q
    const arma::mat deta = d.X * dB + dM * C.t() + M * dC.t()
                         + 0.5 * (dS2 * C2.t() + S2 * dC2.t());   // n×p
    const arma::mat dA   = A % deta;
    const arma::mat AmY  = A - d.Y;

    dgB  = Xw.t() * dA;                                                   // Xwᵀ dA
    dgM  = dA * C + AmY * dC + dM;                       dgM.each_col()  %= w;
    dgPS = 0.5 * (dS2 % (1.0 + A * C2) + S2 % (dA * C2 + A * dC2));
                                                        dgPS.each_col() %= w;
    arma::mat AmYw = AmY; AmYw.each_col() %= w;
    arma::mat dAw  = dA;  dAw.each_col()  %= w;
    arma::mat S2w  = S2;  S2w.each_col()  %= w;
    arma::mat dS2w = dS2; dS2w.each_col() %= w;
    dgC  = dAw.t() * M + AmYw.t() * dM
         + (dA.t() * S2w + A.t() * dS2w) % C
         + (A.t() * S2w) % dC;
}

// Diagonal of the outer-block Hessian L_θθ (analytic, strictly positive) — Jacobi
// preconditioner for the CG. Ignores the Schur correction (H_red ≤ L_θθ) but captures
// the dominant B-vs-C curvature scale that otherwise makes the CG ill-conditioned.
//   diag(L_BB)_{ab} = Σ_i w_i X_ia² A_ib
//   diag(L_CC)_{jk} = Σ_i w_i A_ij [(M_ik + S²_ik C_jk)² + S²_ik]
inline void precond_diag(const PlnData & d, const arma::mat & XX,   // XX = X⊙X
                         const arma::mat & C, const arma::mat & C2, const arma::mat & M,
                         const arma::mat & A, const arma::mat & S2,
                         arma::mat & pB, arma::mat & pC) {
    arma::mat Aw = A; Aw.each_col() %= d.w;
    pB = XX.t() * Aw;                                              // d×p
    pC = Aw.t() * (M % M) + 2.0 * (C % (Aw.t() * (M % S2)))
       + C2 % (Aw.t() * (S2 % S2)) + Aw.t() * S2;                 // p×q
    pB = arma::clamp(pB, 1e-8, arma::datum::inf);
    pC = arma::clamp(pC, 1e-8, arma::datum::inf);
}

} // namespace trnewton
