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

template<typename Traits>
Rcpp::List newton_optimize_impl(
    const arma::mat & Y, const arma::mat & X, const arma::mat & O, const arma::vec & w,
    arma::mat B, arma::mat M, arma::mat S,
    typename Traits::State state,
    int maxiter, double ftol, int max_em, double em_tol,
    int block_newton_thresh = 30
) {
    const int n          = Y.n_rows;
    const arma::uword p  = Y.n_cols;
    const double w_bar   = arma::accu(w);
    const double c1      = 1e-4;
    const arma::mat ones_row = arma::ones(n, 1);

    // Precompute X'WX and P_X once: P_X = (X'WX)^{-1}X'W (d×n) for live B = P_X * M_full
    const arma::mat Xw   = X.each_col() % w;    // n×d
    const arma::mat XtWX = X.t() * Xw;           // d×d, symmetric PD
    // When d=0 (no covariates), X'WX is 0×0: skip solve to avoid spurious singularity warning
    const arma::mat P_X  = (X.n_cols > 0) ? arma::solve(XtWX, Xw.t()) : arma::mat(0, n);

    arma::mat psi = arma::log(S % S);
    arma::mat S2  = arma::exp(psi);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iter   = 0;
    int last_status  = 5;

    // Inner loop: B profiled at each Newton step via envelope theorem.
    auto inner_loop = [&]() {
        double obj_prev = arma::datum::inf;
        // S2 is current (fixed_point_psi updated it in the previous round, or initialized above).
        // Sync B, Z, A once per EM round; they are kept current at the end of each Newton step
        // so subsequent iterations reuse them directly without recomputation.
        B            = P_X * M;
        arma::mat Z  = O + M;
        arma::mat A  = arma::exp(Z + 0.5 * S2);

        for (int it = 0; it < maxiter; it++) {
            // S2, B, Z, A are all current here.
            const arma::mat XB    = X * B;         // frozen during the Armijo line search below
            const arma::mat M_res = M - XB;

            // Newton step for M_full using Trait-specific preconditioner.
            // DenseOmegaImpl: shared p×p preconditioner (diag(wA_mean)+Omega), one Cholesky per call.
            // Diagonal/Spherical: diagonal step (unchanged).
            arma::mat grad_M, step_M;
            Traits::compute_step_M(M_res, state, A, Y, w, ones_row, grad_M, step_M,
                                   static_cast<arma::uword>(block_newton_thresh));

            // Q_step = (I - X*P_X)*step_M: change in M_res per unit alpha (correct projected step)
            // With live B: M_res_t = M_res - alpha*Q_step, slope = -grad_M . Q_step
            const arma::mat Q_step  = step_M - X * (P_X * step_M);
            arma::mat MresO  = Traits::times_Omega(M_res, state);
            arma::mat QstepO = Traits::times_Omega(Q_step, state);
            double f0_M    = arma::accu(w.t() * (A - Y % Z)) + Traits::penalty_M(MresO, M_res, w);
            double slope_M = -arma::accu(grad_M % Q_step);
            // Fall back if step is not a descent direction (degenerate case)
            if (slope_M >= 0) slope_M = -arma::accu(grad_M % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat MresOt = MresO  - alpha_M * QstepO;
                arma::mat MresT  = M_res  - alpha_M * Q_step;
                arma::mat Zt     = Z      - alpha_M * step_M;  // = O + Mt
                arma::mat At     = arma::exp(Zt + 0.5 * S2);
                if (arma::accu(w.t() * (At - Y % Zt)) + Traits::penalty_M(MresOt, MresT, w)
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            B  = P_X * M;
            Z  = O + M;

            // S step: exact fixed-point ψ = −log(A + diag_Ω)
            fixed_point_psi(psi, S2, Z, A, Traits::cov_diag(state, ones_row));

            // A with new S2 — needed for obj computation and for next iteration
            A = arma::exp(Z + 0.5 * S2);
            const arma::mat M_res_new = M - X * B;
            double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                       + Traits::objective_cov(M_res_new, S2, state, w);
            objective_vec.push_back(obj);
            total_iter++;

            if (it > 0 && converged(obj, obj_prev, ftol)) { last_status = 3; break; }
            obj_prev = obj;
        }
    };

    if (Traits::has_em) {
        for (int em = 0; em < max_em; em++) {
            inner_loop();

            // S2 is current: fixed_point_psi updated it inside inner_loop.
            // B is already at its optimum (live-updated in inner_loop); only Omega needs updating.
            // Recompute M_res for the Sigma/Omega M-step.
            arma::mat M_res = M - X * B;
            Traits::mstep(state, M_res, S2, w, w_bar, p);

            arma::mat Z = O + M;
            arma::mat A = arma::exp(Z + 0.5 * S2);
            double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * psi))
                        + Traits::elbo_cov(state, w_bar, p);
            if (em > 0 && converged(elbo, elbo_prev, em_tol)) { last_status = 3; break; }
            elbo_prev = elbo;
        }
    } else {
        // Fixed covariance: Omega is fixed; B is live-updated in inner_loop.
        double elbo_prev_fix = -arma::datum::inf;
        for (int em = 0; em < max_em; em++) {
            inner_loop();
            // S2 is current: fixed_point_psi updated it inside inner_loop.
            // B already optimal after inner_loop; recompute M_res for ELBO check.
            arma::mat M_res = M - X * B;
            arma::mat Z     = O + M;
            arma::mat A     = arma::exp(Z + 0.5 * S2);
            double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * psi))
                        - Traits::objective_cov(M_res, S2, state, w);
            if (em > 0 && converged(elbo, elbo_prev_fix, em_tol)) { last_status = 3; break; }
            elbo_prev_fix = elbo;
        }
    }

    S2 = arma::exp(psi);
    S  = arma::exp(0.5 * psi);
    arma::mat M_res = M - X * B;
    arma::mat Z = O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    // final_loglik uses M_res for the KL terms (same convention as newton_impl.h)
    arma::vec loglik = Traits::final_loglik(Y, Z, A, M_res, psi, state);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    Rcpp::List cov_out = Traits::output_cov(M_res, S2, w, w_bar, state);
    return Rcpp::List::create(
        Rcpp::Named("B",     B               ),
        Rcpp::Named("M",     M               ),   // M_full
        Rcpp::Named("S",     S               ),
        Rcpp::Named("Z",     Z               ),
        Rcpp::Named("A",     A               ),
        Rcpp::Named("Sigma", cov_out["Sigma"]),
        Rcpp::Named("Omega", cov_out["Omega"]),
        Rcpp::Named("Ji",    Ji              ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status   ),
            Rcpp::Named("backend",    "newton"  ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", total_iter    )
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
    int maxiter, double ftol,
    int block_newton_thresh = 30
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

    // Z and A are kept current at the end of each iteration; compute them once here.
    arma::mat Z = O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    for (int it = 0; it < maxiter; it++) {
        // S2, Z, A are current.
        arma::mat M_res = M - XB;          // M_res for KL terms (B is fixed)

        // ---- Newton step for M (gradient and preconditioned step via Trait method) ----
        arma::mat grad_M, step_M;
        Traits::compute_step_M(M_res, state, A, Y, w, ones_row, grad_M, step_M,
                               static_cast<arma::uword>(block_newton_thresh));
        arma::mat MO   = Traits::times_Omega(M_res, state);
        arma::mat dMO  = Traits::times_Omega(step_M, state);
        double f0_M    = arma::accu(w.t() * (A - Y % Z)) + Traits::penalty_M(MO, M_res, w);
        double slope_M = -arma::accu(grad_M % step_M);
        double alpha_M = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            arma::mat MresT = M_res - alpha_M * step_M;
            arma::mat MOt   = MO   - alpha_M * dMO;
            arma::mat Zt    = Z    - alpha_M * step_M;
            arma::mat At    = arma::exp(Zt + 0.5 * S2);
            if (arma::accu(w.t() * (At - Y % Zt)) + Traits::penalty_M(MOt, MresT, w)
                <= f0_M + c1 * alpha_M * slope_M) break;
            alpha_M *= 0.5;
        }
        M -= alpha_M * step_M;
        Z  = O + M;

        // ---- S step: ψ = −log(A + ω²) — exact minimiser for fixed A ----
        fixed_point_psi(psi, S2, Z, A, Traits::cov_diag(state, ones_row));

        // ---- Objective for convergence ----
        A = arma::exp(Z + 0.5 * S2);
        arma::mat M_res_new = M - XB;
        double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                   + Traits::objective_cov(M_res_new, S2, state, w);
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
    }

    // ---- Final output ----
    // S2, Z, A are current from the last iteration (fixed_point_psi + exp update).
    S  = arma::exp(0.5 * psi);
    const arma::mat M_res = M - XB;      // for loglik KL terms
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
