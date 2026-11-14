#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_zipln.h"

// ─────────────────────────────────────────────────────────────────────────────
// Builtin Newton VE-step for ZIPLN: joint (M, ψ = log S², R), Omega fixed.
//
// Same as the inner loop of builtin_optimize_pln_impl (B profiled, joint 2x2
// Newton step), plus:
//   - R set to its conditional optimum R* = σ(A + logit(Pi)) where Y = 0, else 0,
//     at the top of every iteration, and frozen during the Newton step and its
//     line search (∂f/∂R = 0 at R*);
//   - A replaced by A_eff = (1-R) ⊙ A (see covariance_zipln.h);
//   - unit weights (ZIPLN has no observation weights).
// B is re-profiled at every Newton step, which the joint (M, ψ, R) step needs
// to track the optimum.
//
// cg_max > 1: truncated Newton step (Newton-CG). The cell-wise 2×2 step ignores
// the coupling of the cells through Omega (off its diagonal) and through the
// profiling of B (projection Q = I - X P_X). It is slow on the zero counts,
// where A_eff ≈ 0 leaves the diagonal of the Hessian far from dominant. The
// Newton system in (M, ψ), coupling included, is solved by at most cg_max
// conjugate gradient iterations, preconditioned by the 2×2 blocks (the first
// iterate is the cell-wise step, up to its length), stopped at a relative
// residual of cg_tol. Each iteration costs one product by Omega; the step and
// its product by Omega are combinations of the CG directions and of their
// products, so the line search stays free of products. Unit weights assumed
// (the Hessian of the penalty is then v -> Q v Omega).
template <typename Traits>
Rcpp::List builtin_vestep_zipln_impl(
    const PlnData & d,
    const arma::mat & init_M, const arma::mat & init_S2,
    const arma::mat & Pi, const arma::mat & B,
    const typename Traits::State & state,
    int maxiter, double ftol_rel, int cg_max = 1, double cg_tol = 0.1
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

    // V -> Q V, projection on null(X') (identity without covariates)
    auto project = [&](const arma::mat & V) -> arma::mat {
        return do_profile ? arma::mat(V - d.X * (P_X * V)) : V;
    };

    // Truncated Newton step by preconditioned conjugate gradients on
    // H (step_M, step_psi) = (grad_M, grad_psi), with
    //   H_mm v = A_eff ⊙ v + Q v Omega,  H_mp = h_mp (cell-wise),  H_pp = h_pp (cell-wise)
    // (Q (Q v Omega) = Q v Omega since Q v Omega is in null(X')).
    // Also returns Q_step and Q_step Omega, accumulated along the CG iterations.
    auto newton_cg_step = [&](const arma::mat & A_eff, const arma::mat & S2,
                              arma::mat & grad_M, arma::mat & grad_psi, arma::mat & MresO,
                              arma::mat & step_M, arma::mat & step_psi,
                              arma::mat & Q_step, arma::mat & QstepO) {
        arma::mat h_mm, h_mp, h_pp, det;
        CovTraitsBase<Traits>::joint_blocks(M_res, state, A_eff, S2, d.Y, d.w,
                                            grad_M, grad_psi, MresO, h_mm, h_mp, h_pp, det);
        step_M.zeros(arma::size(grad_M)); step_psi.zeros(arma::size(grad_M));
        Q_step.zeros(arma::size(grad_M)); QstepO.zeros(arma::size(grad_M));
        arma::mat r_M = grad_M, r_psi = grad_psi, z_M, z_psi;
        CovTraitsBase<Traits>::solve_blocks(h_mm, h_mp, h_pp, det, r_M, r_psi, z_M, z_psi);
        arma::mat s_M = z_M, s_psi = z_psi;
        double rz = arma::accu(r_M % z_M) + arma::accu(r_psi % z_psi);
        const double tol = cg_tol * std::sqrt(arma::accu(arma::square(grad_M)) + arma::accu(arma::square(grad_psi)));
        for (int k = 0; k < cg_max; k++) {
            const arma::mat Qs  = project(s_M);
            const arma::mat QsO = Traits::times_Omega(Qs, state);
            const arma::mat H_M   = A_eff % s_M + QsO + h_mp % s_psi;
            const arma::mat H_psi = h_mp % s_M + h_pp % s_psi;
            const double curv = arma::accu(s_M % H_M) + arma::accu(s_psi % H_psi);
            if (curv <= 0.0) {
                // negative curvature: keep the iterate, or the preconditioned gradient
                if (k == 0) { step_M = s_M; step_psi = s_psi; Q_step = Qs; QstepO = QsO; }
                break;
            }
            const double a = rz / curv;
            step_M += a * s_M; step_psi += a * s_psi; Q_step += a * Qs; QstepO += a * QsO;
            r_M -= a * H_M; r_psi -= a * H_psi;
            if (std::sqrt(arma::accu(arma::square(r_M)) + arma::accu(arma::square(r_psi))) <= tol) break;
            CovTraitsBase<Traits>::solve_blocks(h_mm, h_mp, h_pp, det, r_M, r_psi, z_M, z_psi);
            const double rz_new = arma::accu(r_M % z_M) + arma::accu(r_psi % z_psi);
            s_M   = z_M   + (rz_new / rz) * s_M;
            s_psi = z_psi + (rz_new / rz) * s_psi;
            rz = rz_new;
        }
    };

    // The objective omits the constant sum log(Y!), which dominates it on large counts (e.g.
    // 5.6e7 against a bound of 2e5 on microcosm) and makes a relative tolerance far looser
    // than it reads. The Newton-CG step tests convergence on the objective shifted back by
    // this constant (the cell-wise step keeps its historical test).
    const double f_shift = (cg_max > 1) ? arma::accu(arma::lgamma(d.Y + 1.0)) : 0.0;

    arma::mat R(arma::size(A), arma::fill::zeros);
    std::vector<double> objective_vec;
    int iter = 0;

    for (iter = 0; iter < maxiter; iter++) {
        // Exact conditional optimum of R, frozen for this Newton step + line search.
        R = zipln_update_R(A, ctx);
        const arma::mat one_m_R = 1.0 - R;
        const arma::mat A_eff   = one_m_R % A;

        arma::mat grad_M, step_M, grad_psi, step_psi, MresO, Q_step, QstepO;
        double slope;
        if (cg_max <= 1) {
            CovTraitsBase<Traits>::compute_joint_step_MS(M_res, state, A_eff, S2, d.Y, d.w,
                                                           grad_M, step_M, grad_psi, step_psi, MresO);
            Q_step = project(step_M);
            QstepO = Traits::times_Omega(Q_step, state);
            slope  = -arma::accu(grad_M % Q_step) - arma::accu(grad_psi % step_psi);
        } else {
            newton_cg_step(A_eff, S2, grad_M, grad_psi, MresO, step_M, step_psi, Q_step, QstepO);
            // exact directional derivative: Q M Omega = M_res Omega is in null(X')
            slope  = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
        }
        const double f0 = CovTraitsBase<Traits>::objective(M_res, Z, S2, psi, MresO, state, d.Y, d.w, A_eff);
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
        if (converged(obj + f_shift, f0 + f_shift, ftol_rel)) break;
    }

    return make_zipln_vestep_result(M, S2, R, 3, "newton", objective_vec, iter);
}

// Extracts params/config and calls builtin_vestep_zipln_impl (one instantiation
// per covariance structure in wrappers_builtin_optim_zipln.cpp)
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
    // Newton-CG step (see builtin_vestep_zipln_impl); cg_max = 1 gives the cell-wise step
    const int    cg_max = config.containsElementNamed("cg_max") ? Rcpp::as<int>(config["cg_max"])    : 1;
    const double cg_tol = config.containsElementNamed("cg_tol") ? Rcpp::as<double>(config["cg_tol"]) : 0.1;
    typename Traits::State state(Omega);
    return builtin_vestep_zipln_impl<Traits>(d, init_M, init_S2, Pi, B, state, cfg.maxiter, cfg.ftol, cg_max, cg_tol);
}
