#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "CovarianceTraits.h"

// Generic coordinate-Newton EM optimizer for PLN covariance variants.
// Traits encodes the variant-specific math (M grad/hess, M-step, ELBO, loglik).
// has_em=true → double EM+inner loop; has_em=false → single inner loop (fixed cov).
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

    const arma::mat Xw   = X.each_col() % w;
    arma::mat Xw2        = X % X; Xw2.each_col() %= w;
    const arma::mat ones_row = arma::ones(n, 1);

    arma::mat psi = arma::log(S % S);   // ψ = log(S²), work in logS² space throughout
    arma::mat S2  = arma::exp(psi);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iter   = 0;
    int last_status  = 5;

    auto inner_loop = [&]() {
        double obj_prev = arma::datum::inf;
        // S2 is current (fixed_point_psi updated it in previous round, or initialized above).
        // Z and A are kept current at the end of each iteration; compute them once here.
        arma::mat Z = O + X * B + M;
        arma::mat A = arma::exp(Z + 0.5 * S2);

        for (int it = 0; it < maxiter; it++) {
            // S2, Z, A are current; newton_step_B will update B, Z, A.
            newton_step_B(Xw, Xw2, X, Y, O, w, M, S2, B, Z, A);

            arma::mat grad_M, hess_M;
            Traits::grad_hess_M(M, state, A, Y, w, ones_row, grad_M, hess_M);
            hess_M.clamp(1e-10, arma::datum::inf);
            arma::mat step_M = grad_M / hess_M;
            // Precompute MO and dMO once — avoids O(n*p²) multiply at each Armijo backtrack
            arma::mat MO   = Traits::times_Omega(M, state);
            arma::mat dMO  = Traits::times_Omega(step_M, state);
            double f0_M    = arma::accu(w.t() * (A - Y % Z)) + Traits::penalty_M(MO, M, w);
            double slope_M = -arma::accu(grad_M % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt  = M  - alpha_M * step_M;
                arma::mat MOt = MO - alpha_M * dMO;    // linear update, no matrix multiply
                arma::mat Zt  = Z  - alpha_M * step_M;
                arma::mat At  = arma::exp(Zt + 0.5 * S2);
                if (arma::accu(w.t() * (At - Y % Zt)) + Traits::penalty_M(MOt, Mt, w)
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + X * B + M;

            // S step: ψ = −log(A + ω²)  — exact minimiser for fixed A
            fixed_point_psi(psi, S2, Z, A, Traits::cov_diag(state, ones_row));

            A = arma::exp(Z + 0.5 * S2);
            // KL entropy: S² − log(S²) − 1 = exp(ψ) − ψ − 1
            double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                       + Traits::objective_cov(M, S2, state, w);
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
            Traits::mstep(state, M, S2, w, w_bar, p);

            arma::mat Z = O + X * B + M;
            arma::mat A = arma::exp(Z + 0.5 * S2);
            double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * psi))
                        + Traits::elbo_cov(state, w_bar, p);
            if (em > 0 && converged(elbo, elbo_prev, em_tol)) { last_status = 3; break; }
            elbo_prev = elbo;
        }
    } else {
        inner_loop();
    }

    S2 = arma::exp(psi);
    S  = arma::exp(0.5 * psi);
    arma::mat Z = O + X * B + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec loglik = Traits::final_loglik(Y, Z, A, M, psi, state);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    Rcpp::List cov_out = Traits::output_cov(M, S2, w, w_bar, state);
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
            Rcpp::Named("status",     last_status   ),
            Rcpp::Named("backend",    "newton"      ),
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

    // Z and A are kept current at the end of each iteration; compute them once here.
    arma::mat Z = O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    for (int it = 0; it < maxiter; it++) {
        // S2, Z, A are current.
        arma::mat M_res = M - XB;          // M_res for KL terms (B is fixed)

        // ---- Diagonal Newton step for M (gradient == gradient w.r.t. M_res, B fixed) ----
        arma::mat grad_M, hess_M;
        Traits::grad_hess_M(M_res, state, A, Y, w, ones_row, grad_M, hess_M);
        hess_M.clamp(1e-10, arma::datum::inf);
        arma::mat step_M = grad_M / hess_M;
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
