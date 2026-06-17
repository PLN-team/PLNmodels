#pragma once
#include <RcppArmadillo.h>
#include "nlopt_wrapper.h"
#include "packing.h"
#include "covariance_zipln.h"

// ─────────────────────────────────────────────────────────────────────────────
// nlopt VE-step for ZIPLN: joint (M, ψ = log S²), Omega and B fixed.
//
// R is fixed at its exact conditional optimum R* = σ(A₀ + logit(Pi)) [Y=0],
// computed once from the initial (M, S²) before the nlopt solve — required
// for a consistent (objective, gradient) pair across nlopt's line search.
// The final R is recomputed from the optimised (M, S²) before returning.
//
// B is kept fixed at the value from the preceding M-step (no re-profiling
// inside the nlopt solve, unlike the Newton path: at the start of each VE
// step B is already at its optimum and stays close throughout the small
// nlopt steps, so dynamic profiling here was tested and found to add cost
// without improving the loglik).
//
// Mirrors nlopt_vestep_impl (nlopt_optim_pln.h): same (data, params, config)
// signature, State built from params["Omega"] inside the function.
template <typename Traits>
Rcpp::List ve_step_zipln_nlopt_impl(
    const Rcpp::List & data,    // List(Y, X, O, w)
    const Rcpp::List & params,  // List(M, S2, Pi, B, Omega) — Pi, B, Omega fixed
    const Rcpp::List & config
) {
    const PlnData d(data);
    const auto init_M  = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S2 = Rcpp::as<arma::mat>(params["S2"]);
    const auto Pi      = Rcpp::as<arma::mat>(params["Pi"]);
    const auto B       = Rcpp::as<arma::mat>(params["B"]);
    const auto Omega   = Rcpp::as<arma::mat>(params["Omega"]);
    const typename Traits::State state(Omega);

    const arma::mat psi_init = arma::log(init_S2);
    const auto metadata = tuple_metadata(init_M, psi_init);
    enum { M_ID, PSI_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>  (parameters.data()) = init_M;
    metadata.map<PSI_ID>(parameters.data()) = psi_init;

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB = d.X * B;
    const ZiplnRContext ctx(Pi, d.Y);

    const arma::mat A0      = arma::exp(d.O + init_M + 0.5 * init_S2);
    const arma::mat R0      = zipln_update_R(A0, ctx);
    const arma::mat one_m_R = 1.0 - R0;

    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>  (par);
        const arma::mat psi   = metadata.map<PSI_ID>(par);
        const arma::mat S2    = arma::exp(psi);
        const arma::mat M_res = M - XB;
        const arma::mat Z     = d.O + M;
        const arma::mat A_eff = one_m_R % arma::exp(Z + 0.5 * S2);

        arma::mat gM, gPS;
        const double obj = zipln_vestep_obj_grad<Traits>(M_res, Z, S2, psi, A_eff, state, d.Y, d.w, gM, gPS);
        metadata.map<M_ID>  (grad) = gM;
        metadata.map<PSI_ID>(grad) = gPS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    const arma::mat M    = metadata.copy<M_ID>  (parameters.data());
    const arma::mat psi  = metadata.copy<PSI_ID>(parameters.data());
    const arma::mat S2   = arma::exp(psi);
    const arma::mat A    = arma::exp(d.O + M + 0.5 * S2);
    const arma::mat Rfin = zipln_update_R(A, ctx);

    return make_zipln_vestep_result(M, S2, Rfin,
                                     static_cast<int>(result.status), "nlopt", objective_vec, result.nb_iterations);
}
