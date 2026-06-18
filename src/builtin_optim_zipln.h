#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_zipln.h"

// ─────────────────────────────────────────────────────────────────────────────
// Builtin Newton VE-step for ZIPLN: joint (M, ψ = log S², R), Omega fixed.
//
// Structurally identical to builtin_optimize_pln_impl's inner loop (same
// B-profiling via null(X') projection, same joint 2x2 Newton step from
// CovTraitsBase::compute_joint_step_MS), plus:
//   - R updated to its exact conditional optimum at the top of every outer
//     iteration: R* = σ(A + logit(Pi)) where Y = 0, else 0. Frozen during
//     the Newton step + Armijo line search (∂f/∂R = 0 at R*, so updating R
//     does not invalidate the (M, ψ) gradient/step for that iteration).
//   - A replaced by A_eff = (1-R) ⊙ A everywhere it enters
//     compute_joint_step_MS/objective (see covariance_zipln.h for why this
//     substitution reproduces the ZIPLN formulas exactly).
//   - w = d.w: ZIPLN carries no per-row weights of its own, but threads the
//     PlnData convention through (the R caller passes ones(n)).
//
// B is taken from the preceding M-step and re-profiled at every Newton step
// here (unlike PLN's builtin_vestep_pln_impl, where B/Omega are truly fixed):
// dynamic B-profiling was tested and found necessary for ZIPLN's joint
// (M, ψ, R) VE step to track the moving optimum.
template <typename Traits>
Rcpp::List builtin_vestep_zipln_impl(
    const PlnData & d,
    const arma::mat & init_M, const arma::mat & init_S2,
    const arma::mat & Pi, const arma::mat & B,
    const typename Traits::State & state,
    int maxiter, double ftol_rel
) {
    constexpr double c1 = 1e-4;

    const bool      do_profile = (d.X.n_cols > 0);
    const arma::mat P_X = do_profile
        ? arma::solve(d.X.t() * d.X, d.X.t(), arma::solve_opts::likely_sympd)
        : arma::mat(0, d.Y.n_rows);

    const ZiplnRContext ctx(Pi, d.Y);

    arma::mat M     = init_M;            // M_full
    arma::mat M_res = init_M - d.X * B;  // in null(X') by construction (B = P_X * init_M)
    arma::mat psi   = arma::log(init_S2);
    arma::mat S2    = arma::exp(psi);
    arma::mat Z     = d.O + M;
    arma::mat A     = arma::trunc_exp(Z + 0.5 * S2);

    arma::mat R(arma::size(A), arma::fill::zeros);
    std::vector<double> objective_vec;
    int iter = 0;

    for (iter = 0; iter < maxiter; iter++) {
        // Exact conditional optimum of R, frozen for this Newton step + line search.
        R = zipln_update_R(A, ctx);
        const arma::mat one_m_R = 1.0 - R;
        const arma::mat A_eff   = one_m_R % A;

        arma::mat grad_M, step_M, grad_psi, step_psi, MresO;
        CovTraitsBase<Traits>::compute_joint_step_MS(M_res, state, A_eff, S2, d.Y, d.w,
                                                       grad_M, step_M, grad_psi, step_psi, MresO);
        const double f0 = CovTraitsBase<Traits>::objective(M_res, Z, S2, psi, MresO, state, d.Y, d.w, A_eff);

        const arma::mat Q_step = do_profile ? step_M - d.X * (P_X * step_M) : step_M;
        const arma::mat QstepO = Traits::times_Omega(Q_step, state);

        double slope = -arma::accu(grad_M % Q_step) - arma::accu(grad_psi % step_psi);
        if (slope >= 0.0) slope = -(arma::accu(arma::square(grad_M)) + arma::accu(arma::square(grad_psi)));

        double alpha = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            const arma::mat MresT  = M_res - alpha * Q_step;
            const arma::mat MresOt = MresO - alpha * QstepO;
            const arma::mat psit   = psi   - alpha * step_psi;
            const arma::mat S2t    = arma::exp(psit);
            const arma::mat Zt     = Z     - alpha * step_M;
            const arma::mat AeffT  = one_m_R % arma::trunc_exp(Zt + 0.5 * S2t);
            if (CovTraitsBase<Traits>::objective(MresT, Zt, S2t, psit, MresOt, state, d.Y, d.w, AeffT) <= f0 + c1 * alpha * slope) break;
            alpha *= 0.5;
        }

        M_res -= alpha * Q_step;
        const arma::mat MO = MresO - alpha * QstepO;
        M     -= alpha * step_M;
        psi   -= alpha * step_psi;
        S2     = arma::exp(psi);
        Z      = d.O + M;
        A      = arma::trunc_exp(Z + 0.5 * S2);

        const double obj = CovTraitsBase<Traits>::objective(M_res, Z, S2, psi, MO, state, d.Y, d.w, one_m_R % A);
        objective_vec.push_back(obj);
        if (converged(obj, f0, ftol_rel)) break;
    }

    return make_zipln_vestep_result(M, S2, R, 3, "newton", objective_vec, iter);
}

// Thin wrapper: extracts params/config (mirroring builtin_vestep_pln_impl's wrappers)
// and delegates to builtin_vestep_zipln_impl. One instantiation per covariance structure
// in wrappers_builtin_optim_zipln.cpp — each export (builtin_optimize_vestep_zipln_<type>)
// is a one-line call to this template.
template <typename Traits>
Rcpp::List builtin_optimize_vestep_zipln_wrapper(
    const Rcpp::List & data, const Rcpp::List & params, const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const NewtonConfig cfg(config);
    typename Traits::State state(Omega);
    return builtin_vestep_zipln_impl<Traits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol);
}
