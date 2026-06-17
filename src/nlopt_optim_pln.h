#pragma once
#include <RcppArmadillo.h>
#include "nlopt_wrapper.h"
#include "packing.h"
#include "covariance_pln.h"            // CovTraits + ki() (via utils.h)

// ─── Generic EM optimizer: inner nlopt with Omega fixed, outer M-step on Omega ─
// CovTraits must provide:
// State::update(M_res, S2, w, w_bar) [mstep], vestep_obj_grad,
// elbo_cov, final_loglik, output_cov, has_em.
// FullCovTraits  (has_em=true):  loops until ELBO converges, M-step updates Omega.
// FixedCovTraits (has_em=false): single inner nlopt pass, Omega fixed externally.
// The caller constructs the initial State (4-arg ctor from data for Full,
// 1-arg ctor from external Omega for Fixed) since that's the one part that
// genuinely differs between the two use cases.
template <typename CovTraits>
Rcpp::List nlopt_optimize_em_impl(
    const Rcpp::List & data,    // List(Y, X, O, w)
    const arma::mat & init_M, const arma::mat & init_S2,
    typename CovTraits::State state,
    const Rcpp::List & config
) {
    const PlnData d(data);

    const arma::uword p = d.Y.n_cols;
    const double w_bar  = arma::accu(d.w);
    const NewtonConfig cfg(config);

    const arma::mat Xw  = d.X.each_col() % d.w;
    const arma::mat P_X = (d.X.n_cols > 0) ? arma::solve(d.X.t() * Xw, Xw.t()) : arma::mat(0, d.Y.n_rows);

    const auto metadata = tuple_metadata(init_M, init_S2);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S2);

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iterations = 0;
    int last_status = 0;

    for (int em_iter = 0; em_iter < std::max(1, cfg.max_em); em_iter++) {
        auto optimizer = new_nlopt_optimizer(config, parameters.size());
        objective_vec.reserve(objective_vec.size() + nlopt_get_maxeval(optimizer.get()));

        // Inner E-step: M_full is the NLOPT parameter; B profiled at each eval
        // (envelope theorem); Omega fixed at the current state for the duration
        // of this inner optimization.
        auto objective_and_grad = [&](const double * par, double * grad) -> double {
            const arma::mat M_full = metadata.map<M_ID>(par);
            const arma::mat logS2  = metadata.map<S_ID>(par);
            const arma::mat S2     = arma::exp(logS2);
            const arma::mat B      = P_X * M_full;
            const arma::mat M_res  = M_full - d.X * B;
            arma::mat gM, gS;
            const double obj = CovTraits::vestep_obj_grad(M_res, d.O + M_full, S2, logS2, state, d.Y, d.w, gM, gS);
            metadata.map<M_ID>(grad) = gM;
            metadata.map<S_ID>(grad) = gS;
            objective_vec.push_back(obj);
            return obj;
        };

        OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
        total_iterations += result.nb_iterations;
        last_status = static_cast<int>(result.status);

        if constexpr (CovTraits::has_em) {
            arma::mat M_full = metadata.copy<M_ID>(parameters.data());
            arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
            arma::mat S2     = arma::exp(logS2);
            arma::mat B      = P_X * M_full;
            arma::mat M_res  = M_full - d.X * B;
            state.update(M_res, S2, d.w, w_bar);  // M-step

            arma::mat Z = d.O + M_full;
            arma::mat A = arma::exp(Z + 0.5 * S2);
            double elbo = arma::accu(d.w.t() * (d.Y % Z - A + 0.5 * logS2)) + CovTraits::elbo_cov(state, w_bar, p);
            if (em_iter > 0 && converged(elbo, elbo_prev, cfg.em_tol)) break;
            elbo_prev = elbo;
        } else {
            break;  // fixed covariance: single inner pass, no M-step
        }
    }

    arma::mat M      = metadata.copy<M_ID>(parameters.data());  // M_full
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat B      = P_X * M;
    arma::mat M_res  = M - d.X * B;
    arma::mat Z      = d.O + M;
    arma::mat A      = arma::exp(Z + 0.5 * S2);

    arma::vec loglik = CovTraits::final_loglik(d.Y, Z, A, M_res, logS2, state);
    Rcpp::List cov_out = CovTraits::output_cov(M_res, S2, d.w, w_bar, state);
    return make_pln_result(B, M, S2, Z, A, cov_out, loglik, last_status, "nlopt", objective_vec, total_iterations);
}

// ─── Generic joint optimizer with Omega profiled at every eval ────────────────
// CovTraits must provide:
//   - State(M_res, S2, w, w_bar)  — 4-arg constructor that profiles Omega from data
//   - static profile_and_grad(s, M_res, Z, S2, logS2, Y, w, w_bar, p, gM, gPS)
//   - static final_loglik, output_cov, State::update
// Covers DiagonalCovTraits and SphericalCovTraits (no EM loop needed).
// Also covers FullCovTraits for small p (via nlopt_optimize_full_profiled).
template <typename CovTraits>
Rcpp::List nlopt_joint_profiled_impl(
    const Rcpp::List & data,    // List(Y, X, O, w)
    const Rcpp::List & params,  // List(B [unused], M, S2)
    const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M   = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2  = Rcpp::as<arma::mat>(params["S2"]);

    const arma::uword p = d.Y.n_cols;
    const double w_bar  = arma::accu(d.w);

    const arma::mat Xw  = d.X.each_col() % d.w;
    const arma::mat P_X = (d.X.n_cols > 0) ? arma::solve(d.X.t() * Xw, Xw.t()) : arma::mat(0, (arma::uword)d.Y.n_rows);

    // Initialize state by profiling Omega from the initial (M_res, S2)
    const arma::mat B_init     = P_X * init_M;
    const arma::mat M_res_init = init_M - d.X * B_init;
    typename CovTraits::State state(M_res_init, init_S2, d.w, w_bar);

    const auto metadata = tuple_metadata(init_M, init_S2);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S2);

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    // Joint E-step: B and Omega are profiled analytically at every function evaluation.
    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M_full = metadata.map<M_ID>(par);
        const arma::mat logS2  = metadata.map<S_ID>(par);
        const arma::mat S2     = arma::exp(logS2);
        const arma::mat B      = P_X * M_full;
        const arma::mat M_res  = M_full - d.X * B;
        arma::mat gM, gS;
        const double obj = CovTraits::profile_and_grad(state, M_res, d.O + M_full, S2, logS2,
                                                        d.Y, d.w, w_bar, p, gM, gS);
        metadata.map<M_ID>(grad) = gM;
        metadata.map<S_ID>(grad) = gS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M     = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2    = arma::exp(logS2);
    arma::mat B     = P_X * M;
    arma::mat M_res = M - d.X * B;
    state.update(M_res, S2, d.w, w_bar);  // ensure state matches final parameters

    arma::mat Z = d.O + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);

    arma::vec loglik = CovTraits::final_loglik(d.Y, Z, A, M_res, logS2, state);
    Rcpp::List cov_out = CovTraits::output_cov(M_res, S2, d.w, w_bar, state);
    return make_pln_result(B, M, S2, Z, A, cov_out, loglik,
                            static_cast<int>(result.status), "nlopt", objective_vec, result.nb_iterations);
}

// ─── Generic VE-step optimizer (M, ψ=log S² only; B and Ω fixed) ─────────────
// CovTraits must provide:
//   - State(const arma::mat & Omega_mat)  — explicit constructor from dense Omega
//   - static vestep_obj_grad(M_res, Z, S2, logS2, State, Y, w, gM, gPS)
//   - static final_loglik(Y, Z, A, M_res, logS2, State)
template <typename CovTraits>
Rcpp::List nlopt_vestep_impl(
    const Rcpp::List & data,    // List(Y, X, O, w)
    const Rcpp::List & params,  // List(M, S2, B, Omega) — B, Omega fixed
    const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M   = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2  = Rcpp::as<arma::mat>(params["S2"]);
    const auto B        = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega    = Rcpp::as<arma::mat>(params["Omega"]);

    const auto metadata = tuple_metadata(init_M, init_S2);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S2);

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB = d.X * B;
    const typename CovTraits::State s(Omega);

    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(par);
        const arma::mat logS2 = metadata.map<S_ID>(par);
        const arma::mat S2    = arma::exp(logS2);
        const arma::mat M_res = M - XB;
        arma::mat gM, gS;
        const double obj = CovTraits::vestep_obj_grad(M_res, d.O + M, S2, logS2, s, d.Y, d.w, gM, gS);
        metadata.map<M_ID>(grad) = gM;
        metadata.map<S_ID>(grad) = gS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    const arma::mat M     = metadata.copy<M_ID>(parameters.data());
    const arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    const arma::mat S2    = arma::exp(logS2);
    const arma::mat M_res = M - XB;
    const arma::mat Z     = d.O + M;
    const arma::mat A     = arma::exp(Z + 0.5 * S2);

    arma::vec loglik = CovTraits::final_loglik(d.Y, Z, A, M_res, logS2, s);
    return make_vestep_result(M, S2, loglik, static_cast<int>(result.status), "nlopt", objective_vec, result.nb_iterations);
}
