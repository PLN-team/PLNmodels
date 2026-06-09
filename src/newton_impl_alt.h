#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "CovarianceTraits.h"

// Alternative parameterization: M stores M_full (the full variational mean on Z_i,
// not the residual M_res = M_full - XB used in newton_impl.h).
//
// Key difference from the standard implementation:
//   - B is NOT updated in the inner loop; it gets an exact closed-form EM M-step:
//       B = (X'WX)^{-1} X'W M_full
//   - The inner loop updates only M_full and S, with gradient using M_res = M_full - XB.
//   - Sigma/Omega M-step uses M_res = M_full - XB  (same formula as before).
//
// At input/output: M is in the "residual" format (M_res = M_full - XB) to stay
// compatible with the rest of the R package.  The conversion is done internally.

template<typename Traits>
Rcpp::List newton_optimize_alt_impl(
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

    // Precompute X'WX once: used for closed-form B = (X'WX)^{-1} X'W M_full
    const arma::mat Xw   = X.each_col() % w;    // n×d
    const arma::mat XtWX = X.t() * Xw;           // d×d, symmetric PD

    // Convert input M from residual format to M_full = XB + M_res
    M = X * B + M;

    arma::mat psi = arma::log(S % S);
    arma::mat S2  = arma::exp(psi);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iter   = 0;
    int last_status  = 5;

    // Inner loop: B is fixed; only M_full and S are updated.
    auto inner_loop = [&]() {
        const arma::mat XB = X * B;   // constant within this inner loop
        double obj_prev = arma::datum::inf;

        for (int it = 0; it < maxiter; it++) {
            S2 = arma::exp(psi);
            arma::mat M_res = M - XB;
            arma::mat Z = O + M;               // Z = O + M_full
            arma::mat A = arma::exp(Z + 0.5 * S2);

            // Newton step for M_full: gradient/Hessian are functions of M_res
            arma::mat grad_M, hess_M;
            Traits::grad_hess_M(M_res, state, A, Y, w, ones_row, grad_M, hess_M);
            hess_M.clamp(1e-10, arma::datum::inf);
            arma::mat step_M = grad_M / hess_M;

            arma::mat MresO = Traits::times_Omega(M_res, state);
            arma::mat dMO   = Traits::times_Omega(step_M, state);
            double f0_M    = arma::accu(w.t() * (A - Y % Z)) + Traits::penalty_M(MresO, M_res, w);
            double slope_M = -arma::accu(grad_M % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt     = M   - alpha_M * step_M;  // M_full candidate
                arma::mat MresOt = MresO - alpha_M * dMO;
                arma::mat Zt     = Z   - alpha_M * step_M;  // = O + Mt
                arma::mat At     = arma::exp(Zt + 0.5 * S2);
                // penalty uses M_res_t = Mt - XB
                if (arma::accu(w.t() * (At - Y % Zt)) + Traits::penalty_M(MresOt, Mt - XB, w)
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + M;

            // S step: exact fixed-point ψ = −log(A + diag_Ω)
            fixed_point_psi(psi, S2, Z, A, Traits::cov_diag(state, ones_row));

            A = arma::exp(Z + 0.5 * S2);
            M_res = M - XB;
            double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * psi))
                       + Traits::objective_cov(M_res, S2, state, w);
            objective_vec.push_back(obj);
            total_iter++;

            if (it > 0 && converged(obj, obj_prev, ftol)) { last_status = 3; break; }
            obj_prev = obj;
        }
    };

    if (Traits::has_em) {
        for (int em = 0; em < max_em; em++) {
            inner_loop();

            S2 = arma::exp(psi);

            // Closed-form B M-step: B = (X'WX)^{-1} X'W M_full
            B = arma::solve(XtWX, Xw.t() * M);

            // Sigma/Omega M-step: uses M_res = M_full - XB (same as current mstep)
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
        // Fixed covariance: Omega is fixed; iterate (E-step + B M-step) until convergence
        double elbo_prev_fix = -arma::datum::inf;
        for (int em = 0; em < max_em; em++) {
            inner_loop();
            S2 = arma::exp(psi);
            B  = arma::solve(XtWX, Xw.t() * M);
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
        Rcpp::Named("M",     M_res           ),   // return residual format for R compatibility
        Rcpp::Named("S",     S               ),
        Rcpp::Named("Z",     Z               ),
        Rcpp::Named("A",     A               ),
        Rcpp::Named("Sigma", cov_out["Sigma"]),
        Rcpp::Named("Omega", cov_out["Omega"]),
        Rcpp::Named("Ji",    Ji              ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status   ),
            Rcpp::Named("backend",    "newton_alt"  ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", total_iter    )
        ))
    );
}
