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

    arma::mat S2   = S % S;
    arma::mat logS = arma::log(S);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iter   = 0;
    int last_status  = 5;

    auto inner_loop = [&]() {
        double obj_prev = arma::datum::inf;
        for (int it = 0; it < maxiter; it++) {
            S2 = S % S;
            arma::mat Z = O + X * B + M;
            arma::mat A = arma::exp(Z + 0.5 * S2);

            newton_step_B(Xw, Xw2, X, Y, O, w, M, S2, B, Z, A);

            arma::mat grad_M, hess_M;
            Traits::grad_hess_M(M, state, A, Y, w, ones_row, grad_M, hess_M);
            hess_M.clamp(1e-10, arma::datum::inf);
            arma::mat step_M = grad_M / hess_M;
            double f0_M    = arma::accu(w.t() * (A - Y % Z)) + Traits::penalty_M(M, state, w);
            double slope_M = -arma::accu(grad_M % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt = M - alpha_M * step_M;
                arma::mat Zt = Z - alpha_M * step_M;
                arma::mat At = arma::exp(Zt + 0.5 * S2);
                if (arma::accu(w.t() * (At - Y % Zt)) + Traits::penalty_M(Mt, state, w)
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + X * B + M;

            fixed_point_logS(logS, S, S2, Z, A, Traits::cov_diag(state, ones_row));

            A = arma::exp(Z + 0.5 * S2);
            double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * arma::trunc_log(S2)))
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

            S2 = S % S;
            Traits::mstep(state, M, S2, w, w_bar, p);

            arma::mat Z = O + X * B + M;
            arma::mat A = arma::exp(Z + 0.5 * S2);
            double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * arma::trunc_log(S2)))
                        + Traits::elbo_cov(state, w_bar, p);
            if (em > 0 && converged(elbo, elbo_prev, em_tol)) { last_status = 3; break; }
            elbo_prev = elbo;
        }
    } else {
        inner_loop();
    }

    S2 = S % S;
    arma::mat Z = O + X * B + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec loglik = Traits::final_loglik(Y, Z, A, M, S2, state);

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
