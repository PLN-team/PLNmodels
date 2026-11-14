#pragma once
#include <RcppArmadillo.h>
#include "utils.h"
#include "covariance_pln.h"

// M is the full variational mean of Z_i (= X_i*B + M_res), on input and output.
// B is profiled at every Newton step: B = P_X * M = (X'WX)^{-1} X'W M, and the
// gradient w.r.t. M is the one w.r.t. M_res = M - X*B (envelope theorem).
//
// Outer EM loop (max_em) over an inner VE-step optimizing (M, ψ) jointly at fixed
// Omega, then one M-step for Omega. The VE step is a Newton method: cell-wise 2×2
// step (opts.cg_max <= 1), or truncated Newton step by preconditioned conjugate
// gradients (opts.cg_max > 1, see CovTraitsBase::newton_cg_step).
//
// The EM loop converges linearly, slowly when the latent means of small counts are
// poorly determined (rate 0.994 on microcosm). With opts.accelerate, it is accelerated
// by SQUAREM (Varadhan and Roland, 2008, scheme S3, adaptive bound on the step length):
// the map F is one EM iteration (VE step, then M step) on the state (M, ψ), B being
// profiled; two iterations, then one from the extrapolated state, kept if it does
// better. max_em counts the EM iterations (VE steps).
template<typename Traits>
Rcpp::List builtin_optimize_pln_impl(
    const PlnData & d,
    arma::mat B, arma::mat M, arma::mat S2,
    typename Traits::State state,
    int maxiter, double ftol, int max_em, double em_tol, double f_shift = 0.0,
    const SolverOptions & opts = SolverOptions()
) {
    const int n          = d.Y.n_rows;
    const arma::uword p  = d.Y.n_cols;
    const double w_bar   = arma::accu(d.w);
    const double c1      = 1e-4;

    const arma::mat Xw   = d.X.each_col() % d.w;
    const arma::mat XtWX = d.X.t() * Xw;
    const arma::mat P_X  = (d.X.n_cols > 0) ? arma::solve(XtWX, Xw.t()) : arma::mat(0, n);
    auto project = [&](const arma::mat & V) -> arma::mat { return V - d.X * (P_X * V); };

    std::vector<double> objective_vec;
    int last_status  = 5;

    // ── VE step: optimize (M, ψ) at the current Omega/state, to ftol (relative change of
    // the objective, shifted by f_shift). M_res = M - X*B is carried along: B is implicitly
    // re-profiled from the new M at every step (envelope theorem).
    // S2 = exp(psi) is passed along (not recomputed) to keep the cell-wise path bit-exact.
    auto ve_step = [&](arma::mat & M, arma::mat & M_res, arma::mat & psi, arma::mat & S2) {
        arma::mat Z  = d.O + M;
        arma::mat A  = arma::exp(Z + 0.5 * S2);
        double inner_prev = arma::datum::inf;
        for (int it = 0; it < maxiter; it++) {
            arma::mat grad_M, step_M, grad_psi, step_psi, MresO, Q_step, QstepO;
            double slope;
            if (opts.cg_max <= 1) {
                // Joint Newton step: compute_joint_step_MS also returns MO = M_res * Omega
                // so we reuse it for the Armijo penalty — avoids a redundant O(n p²) product.
                Traits::compute_joint_step_MS(M_res, state, A, S2, d.Y, d.w,
                                              grad_M, step_M, grad_psi, step_psi, MresO);
                Q_step = project(step_M);
                QstepO = Traits::times_Omega(Q_step, state);
                slope  = -arma::accu(grad_M % Q_step) - arma::accu(grad_psi % step_psi);
                if (slope >= 0) slope = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
            } else {
                Traits::newton_cg_step(M_res, state, A, S2, d.Y, d.w, project, opts.cg_max, opts.cg_tol,
                                       grad_M, grad_psi, MresO, step_M, step_psi, Q_step, QstepO);
                // exact directional derivative (W M_res Omega is orthogonal to (I - Q) step)
                slope = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
                if (slope >= 0) slope = -(arma::accu(arma::square(grad_M)) + arma::accu(arma::square(grad_psi)));
            }
            const double f0 = Traits::objective(M_res, Z, S2, psi, MresO, state, d.Y, d.w, A);
            double alpha = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                const arma::mat MresT  = M_res - alpha * Q_step;
                const arma::mat MresOt = MresO - alpha * QstepO;
                const arma::mat psit   = psi   - alpha * step_psi;
                const arma::mat S2t    = arma::exp(psit);
                const arma::mat Zt     = Z     - alpha * step_M;
                const arma::mat At     = arma::exp(Zt + 0.5 * S2t);
                if (Traits::objective(MresT, Zt, S2t, psit, MresOt, state, d.Y, d.w, At) <= f0 + c1 * alpha * slope) break;
                alpha *= 0.5;
            }

            // Post-step update: M_res -= alpha*Q_step and MresO -= alpha*QstepO follow from
            // the envelope theorem (B is implicitly re-profiled from the new M), so they're
            // carried forward without ever recomputing B = P_X*M or X*B inside this loop.
            M_res -= alpha * Q_step;
            MresO -= alpha * QstepO;
            M     -= alpha * step_M;
            psi   -= alpha * step_psi;
            S2     = arma::exp(psi);
            Z      = d.O + M;
            A      = arma::exp(Z + 0.5 * S2);

            const double obj = Traits::objective(M_res, Z, S2, psi, MresO, state, d.Y, d.w, A);
            objective_vec.push_back(obj);
            if (it > 0 && converged(obj + f_shift, inner_prev + f_shift, ftol)) break;
            inner_prev = obj;
        }
    };

    // ── Outer ELBO (Omega profiled after the M step), without the constant f_shift
    auto bound = [&](const arma::mat & M, const arma::mat & psi) -> double {
        const arma::mat Z = d.O + M;
        const arma::mat A = arma::exp(Z + 0.5 * arma::exp(psi));   // = A of the last Newton step
        return arma::accu(d.w.t() * (d.Y % Z - A + 0.5 * psi)) + Traits::elbo_cov(state, w_bar, p);
    };

    arma::mat psi   = arma::log(S2);
    // M_res = M - X*B with B = P_X*M; carried across iterations from here on (see above),
    // so B itself only needs to be recomputed once, at the very end for the returned value.
    arma::mat M_res = M - d.X * (P_X * M);

    if (!Traits::has_em) {
        ve_step(M, M_res, psi, S2);
        last_status = 3;   // fixed covariance: inner loop already converged
    } else if (!opts.accelerate) {
        double elbo_prev = -arma::datum::inf;
        for (int em_iter = 0; em_iter < std::max(1, max_em); em_iter++) {
            ve_step(M, M_res, psi, S2);
            Traits::mstep(state, M_res, S2, d.w, w_bar, p);
            const double elbo = bound(M, psi);
            if (em_iter > 0 && converged(elbo - f_shift, elbo_prev - f_shift, em_tol)) { last_status = 3; break; }
            elbo_prev = elbo;
        }
    } else {
        // SQUAREM. The state (Omega) is the M step of the current point on entry of each
        // map evaluation, and of the point returned on exit.
        auto em_map = [&](arma::mat & M, arma::mat & M_res, arma::mat & psi) -> double {
            arma::mat S2_ = arma::exp(psi);
            ve_step(M, M_res, psi, S2_);
            Traits::mstep(state, M_res, S2_, d.w, w_bar, p);
            return bound(M, psi);
        };
        double elbo = em_map(M, M_res, psi);
        int n_em = 1;
        double step_max = 1.0;
        while (n_em < max_em) {
            arma::mat M1 = M, R1 = M_res, P1 = psi;
            const double e1 = em_map(M1, R1, P1);
            arma::mat M2 = M1, R2 = R1, P2 = P1;
            double e_new = em_map(M2, R2, P2);
            n_em += 2;
            const arma::mat rM = M1 - M,  rP = P1 - psi;
            const arma::mat vM = M2 - M1 - rM, vP = P2 - P1 - rP;
            const double nr = arma::accu(arma::square(rM)) + arma::accu(arma::square(rP));
            const double nv = arma::accu(arma::square(vM)) + arma::accu(arma::square(vP));
            const double alpha = (nv > 0.0) ? std::max(1.0, std::min(step_max, std::sqrt(nr / nv))) : 1.0;
            bool keep_2 = true;
            if (alpha > 1.0 && n_em < max_em) {
                arma::mat Mc = M + 2.0 * alpha * rM + alpha * alpha * vM;
                arma::mat Pc = psi + 2.0 * alpha * rP + alpha * alpha * vP;
                arma::mat Rc = project(Mc);
                double ec = -arma::datum::inf;
                if (Pc.is_finite() && Mc.is_finite()) {
                    Traits::mstep(state, Rc, arma::exp(Pc), d.w, w_bar, p);
                    ec = em_map(Mc, Rc, Pc);
                    n_em++;
                }
                if (std::isfinite(ec) && ec >= e_new) {
                    M = Mc; M_res = Rc; psi = Pc; e_new = ec; keep_2 = false;
                    if (alpha == step_max) step_max *= 4.0;
                } else {
                    step_max = std::max(1.0, step_max / 4.0);
                    Traits::mstep(state, R2, arma::exp(P2), d.w, w_bar, p);   // back to the M step of point 2
                }
            } else if (alpha == step_max) step_max *= 4.0;
            if (keep_2) { M = M2; M_res = R2; psi = P2; }
            (void) e1;
            const bool done = converged(e_new - f_shift, elbo - f_shift, em_tol);
            elbo = e_new;
            if (done) { last_status = 3; break; }
        }
    }

    S2 = arma::exp(psi);
    const arma::mat Z = d.O + M;
    const arma::mat A = arma::exp(Z + 0.5 * S2);
    B = P_X * M;  // only needed for the returned value — computed once here

    arma::vec loglik = Traits::final_loglik(d.Y, Z, A, M_res, psi, state);
    Rcpp::List cov_out = Traits::output_cov(M_res, S2, d.w, w_bar, state);
    return make_pln_result(B, M, S2, Z, A, cov_out, loglik, last_status, "newton", objective_vec, (int)objective_vec.size());
}

// ─────────────────────────────────────────────────────────────────────────────
// Generic VE-step Newton optimizer (B and Omega fixed, only M and S updated).
// State must be initialized from a fixed Omega via the explicit constructor.
template<typename Traits>
Rcpp::List builtin_vestep_pln_impl(
    const PlnData & d,
    arma::mat M, arma::mat S2,
    const arma::mat & B, const typename Traits::State & state,
    int maxiter, double ftol, double f_shift = 0.0,
    const SolverOptions & opts = SolverOptions()
) {
    const double c1 = 1e-4;
    const arma::mat XB = d.X * B;

    arma::mat psi = arma::log(S2);

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter  = 0;

    arma::mat Z     = d.O + M;
    arma::mat A     = arma::exp(Z + 0.5 * S2);
    arma::mat M_res = M - XB;  // carried across iterations: XB is fixed, so M_res moves by the same increment as M

    for (int it = 0; it < maxiter; it++) {
        // Joint Newton step: MO = M_res * Omega returned to avoid recomputing it
        // for the Armijo penalty evaluation.
        arma::mat grad_M, step_M, grad_psi, step_psi, MO, dMO;
        if (opts.cg_max <= 1) {
            Traits::compute_joint_step_MS(M_res, state, A, S2, d.Y, d.w,
                                          grad_M, step_M, grad_psi, step_psi, MO);
            dMO = Traits::times_Omega(step_M, state);
        } else {
            // Newton-CG step, B fixed (no projection)
            arma::mat step_copy;
            Traits::newton_cg_step(M_res, state, A, S2, d.Y, d.w, [](const arma::mat & V) { return V; },
                                   opts.cg_max, opts.cg_tol, grad_M, grad_psi, MO, step_M, step_psi, step_copy, dMO);
        }
        double f0    = Traits::objective(M_res, Z, S2, psi, MO, state, d.Y, d.w, A);
        double slope = -arma::accu(grad_M % step_M) - arma::accu(grad_psi % step_psi);
        double alpha = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            const arma::mat MresT = M_res - alpha * step_M;
            const arma::mat MOt   = MO    - alpha * dMO;
            const arma::mat psit  = psi   - alpha * step_psi;
            const arma::mat S2t   = arma::exp(psit);
            const arma::mat Zt    = Z     - alpha * step_M;
            const arma::mat At    = arma::exp(Zt + 0.5 * S2t);
            if (Traits::objective(MresT, Zt, S2t, psit, MOt, state, d.Y, d.w, At) <= f0 + c1 * alpha * slope) break;
            alpha *= 0.5;
        }

        // Post-step update: M_res and MO move by the same alpha*step increments as M
        // (B/XB fixed here, so no re-profiling needed) — carried forward to the next
        // iteration instead of recomputing M - XB from scratch.
        M_res -= alpha * step_M;
        MO    -= alpha * dMO;
        M     -= alpha * step_M;
        psi   -= alpha * step_psi;
        S2     = arma::exp(psi);
        Z      = d.O + M;
        A      = arma::exp(Z + 0.5 * S2);

        double obj = Traits::objective(M_res, Z, S2, psi, MO, state, d.Y, d.w, A);
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj + f_shift, obj_prev + f_shift, ftol)) break;
        obj_prev = obj;
    }

    arma::vec loglik = Traits::final_loglik(d.Y, Z, A, M_res, psi, state);

    return make_vestep_result(M, S2, loglik, 3, "newton", objective_vec, total_iter);
}
