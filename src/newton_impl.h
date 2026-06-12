#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "CovarianceTraits.h"

// M_full parameterization: M is the full variational mean of Z_i (= X_i*B + M_res),
// consistent with the ZIPLN convention.  M_res = M - X*B is computed locally for KL.
//
// B is profiled at every Newton step via the envelope theorem:
//   B = P_X * M = (X'WX)^{-1} X'W M   (closed-form optimum for current M)
//   M_res = M - X*B                     (projection orthogonal to col(X))
// The gradient of J_profiled w.r.t. M equals the gradient w.r.t. M_res (envelope theorem).
//
// Input/output: M is in M_full format throughout.

// Mirror of nlopt_full_cov structure: outer EM loop (max_em) over an inner VE-step that
// optimizes (M, ψ) jointly for fixed Omega until convergence (ftol), then one Omega M-step.
// Joint Newton step per inner iteration: diagonal 2×2 per (i,j) with cross-term H_Mψ.
template<typename Traits>
Rcpp::List newton_optimize_impl(
    const arma::mat & Y, const arma::mat & X, const arma::mat & O, const arma::vec & w,
    arma::mat B, arma::mat M, arma::mat S,
    typename Traits::State state,
    int maxiter, double ftol, int max_em, double em_tol
) {
    const int n          = Y.n_rows;
    const arma::uword p  = Y.n_cols;
    const double w_bar   = arma::accu(w);
    const double c1      = 1e-4;
    const arma::mat ones_row = arma::ones(n, 1);

    const arma::mat Xw   = X.each_col() % w;
    const arma::mat XtWX = X.t() * Xw;
    const arma::mat P_X  = (X.n_cols > 0) ? arma::solve(XtWX, Xw.t()) : arma::mat(0, n);

    arma::mat psi = arma::log(S % S);
    arma::mat S2  = arma::exp(psi);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int last_status  = 5;

    B           = P_X * M;
    arma::mat Z = O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    for (int em_iter = 0; em_iter < std::max(1, max_em); em_iter++) {
        // ── Inner VE-step: optimize (M, ψ) to convergence for the current Omega/state ──
        double inner_prev = arma::datum::inf;
        for (int it = 0; it < maxiter; it++) {
            const arma::mat XB    = X * B;
            const arma::mat M_res = M - XB;

            // Joint Newton step: compute_joint_step_MS also returns MO = M_res * Omega
            // so we reuse it for the Armijo penalty — avoids a redundant O(n p²) product.
            arma::mat grad_M, step_M, grad_psi, step_psi, MresO;
            Traits::compute_joint_step_MS(M_res, state, A, S2, Y, w, ones_row,
                                          grad_M, step_M, grad_psi, step_psi, MresO);
            const arma::mat Q_step  = step_M - X * (P_X * step_M);
            const arma::mat QstepO  = Traits::times_Omega(Q_step, state);
            double f0    = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                         + Traits::penalty_M(MresO, M_res, w) + Traits::penalty_S(S2, state, w);
            double slope = -arma::accu(grad_M % Q_step) - arma::accu(grad_psi % step_psi);
            if (slope >= 0) slope = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
            double alpha = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                const arma::mat MresT  = M_res - alpha * Q_step;
                const arma::mat MresOt = MresO - alpha * QstepO;
                const arma::mat psit   = psi   - alpha * step_psi;
                const arma::mat S2t    = arma::exp(psit);
                const arma::mat Zt     = Z     - alpha * step_M;
                const arma::mat At     = arma::exp(Zt + 0.5 * S2t);
                if (arma::accu(w.t() * (At - Y % Zt - 0.5 * psit))
                    + Traits::penalty_M(MresOt, MresT, w) + Traits::penalty_S(S2t, state, w)
                    <= f0 + c1 * alpha * slope) break;
                alpha *= 0.5;
            }
            M   -= alpha * step_M;
            psi -= alpha * step_psi;
            S2   = arma::exp(psi);
            B    = P_X * M;
            Z    = O + M;
            A    = arma::exp(Z + 0.5 * S2);

            objective_vec.push_back(f0);
            if (it > 0 && converged(f0, inner_prev, ftol)) break;
            inner_prev = f0;
        }

        // ── VM-step: update Omega/Sigma (skipped for fixed covariance) ──
        const arma::mat M_res_cur = M - X * B;
        if (Traits::has_em) {
            Traits::mstep(state, M_res_cur, S2, w, w_bar, p);
        } else {
            last_status = 3; break;   // fixed covariance: inner loop already converged
        }

        // ── Outer ELBO for convergence ──
        double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * psi))
                    + Traits::elbo_cov(state, w_bar, p);
        if (em_iter > 0 && converged(elbo, elbo_prev, em_tol)) { last_status = 3; break; }
        elbo_prev = elbo;
    }

    S2 = arma::exp(psi);
    S  = arma::exp(0.5 * psi);
    arma::mat M_res = M - X * B;
    Z = O + M;
    A = arma::exp(Z + 0.5 * S2);

    arma::vec loglik = Traits::final_loglik(Y, Z, A, M_res, psi, state);
    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    Rcpp::List cov_out = Traits::output_cov(M_res, S2, w, w_bar, state);
    return Rcpp::List::create(
        Rcpp::Named("B",     B               ),
        Rcpp::Named("M",     M               ),
        Rcpp::Named("S",     S               ),
        Rcpp::Named("Z",     Z               ),
        Rcpp::Named("A",     A               ),
        Rcpp::Named("Sigma", cov_out["Sigma"]),
        Rcpp::Named("Omega", cov_out["Omega"]),
        Rcpp::Named("Ji",    Ji              ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status                  ),
            Rcpp::Named("backend",    "newton"                     ),
            Rcpp::Named("objective",  objective_vec                ),
            Rcpp::Named("iterations", (int)objective_vec.size()    )
        ))
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Generic VE-step Newton optimizer (B and Omega fixed, only M and S updated).
// State must be initialized from a fixed Omega via the explicit constructor.
template<typename Traits>
Rcpp::List newton_vestep_impl(
    const arma::mat & Y, const arma::mat & X, const arma::mat & O, const arma::vec & w,
    arma::mat M, arma::mat S,
    const arma::mat & B, const typename Traits::State & state,
    int maxiter, double ftol
) {
    const int n = Y.n_rows;
    const double c1 = 1e-4;
    const arma::mat ones_row = arma::ones(n, 1);
    const arma::mat XB = X * B;

    arma::mat psi = arma::log(S % S);   // ψ = log(S²)
    arma::mat S2  = arma::exp(psi);

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter  = 0;

    arma::mat Z = O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    for (int it = 0; it < maxiter; it++) {
        arma::mat M_res = M - XB;

        // Joint Newton step: MO = M_res * Omega returned to avoid recomputing it
        // for the Armijo penalty evaluation.
        arma::mat grad_M, step_M, grad_psi, step_psi, MO;
        Traits::compute_joint_step_MS(M_res, state, A, S2, Y, w, ones_row,
                                      grad_M, step_M, grad_psi, step_psi, MO);
        const arma::mat dMO = Traits::times_Omega(step_M, state);
        double f0    = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                     + Traits::penalty_M(MO, M_res, w) + Traits::penalty_S(S2, state, w);
        double slope = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
        double alpha = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            const arma::mat MresT = M_res - alpha * step_M;
            const arma::mat MOt   = MO    - alpha * dMO;
            const arma::mat psit  = psi   - alpha * step_psi;
            const arma::mat S2t   = arma::exp(psit);
            const arma::mat Zt    = Z     - alpha * step_M;
            const arma::mat At    = arma::exp(Zt + 0.5 * S2t);
            if (arma::accu(w.t() * (At - Y % Zt - 0.5 * psit))
                + Traits::penalty_M(MOt, MresT, w) + Traits::penalty_S(S2t, state, w)
                <= f0 + c1 * alpha * slope) break;
            alpha *= 0.5;
        }
        M   -= alpha * step_M;
        psi -= alpha * step_psi;
        S2   = arma::exp(psi);
        Z    = O + M;
        A    = arma::exp(Z + 0.5 * S2);

        // Post-update objective: reuse MO_new = MO - alpha*dMO = M_res_new * Omega
        // (exact incremental update, avoids an extra O(n p²) product for full cov).
        const arma::mat MO_new   = MO - alpha * dMO;
        const arma::mat M_res_new = M - XB;
        double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                   + Traits::penalty_M(MO_new, M_res_new, w) + Traits::penalty_S(S2, state, w);
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
    }

    S  = arma::exp(0.5 * psi);
    const arma::mat M_res = M - XB;
    arma::vec loglik = Traits::final_loglik(Y, Z, A, M_res, psi, state);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S")  = S,
        Rcpp::Named("Ji") = Ji,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     3           ),
            Rcpp::Named("backend",    "newton"    ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter  )
        ))
    );
}
