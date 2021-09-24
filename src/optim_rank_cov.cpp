#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List cpp_optimize_rank(
    const Rcpp::List & init_parameters, // List(Theta, B, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_B = Rcpp::as<arma::mat>(init_parameters["B"]);         // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,q)

    const auto metadata = tuple_metadata(init_Theta, init_B, init_M, init_S);
    enum { THETA_ID, B_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<THETA_ID>(packed.data()), per_param_list["Theta"]);
            set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M * B.t();
        arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
        double objective = accu(diagmat(w) * (A - Y % Z)) + 0.5 * accu(diagmat(w) * (M % M + S2 - log(S2) - 1.));

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<B_ID>(grad) = (diagmat(w) * (A - Y)).t() * M + (A.t() * (S2.each_col() % w)) % B;
        metadata.map<M_ID>(grad) = diagmat(w) * ((A - Y) * B + M);
        metadata.map<S_ID>(grad) = diagmat(w) * (S - 1. / S + A * (B % B) % S);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat Sigma = B * (M.t() * (M.each_col() % w) + diagmat(sum(S2.each_col() % w, 0))) * B.t() / accu(w);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M * B.t();
    arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
    arma::mat loglik = arma::sum(Y % Z - A, 1) - 0.5 * sum(M % M + S2 - log(S2) - 1., 1) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", Theta),
        Rcpp::Named("B", B),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// VE rank
// Rank-constrained covariance (for prediction in the PCA space)

// [[Rcpp::export]]
Rcpp::List cpp_optimize_vestep_rank(
        const Rcpp::List & init_parameters, // List(M, S)
        const arma::mat & Y,                // responses (n,p)
        const arma::mat & X,                // covariates (n,d)
        const arma::mat & O,                // offsets (n,p)
        const arma::vec & w,                // weights (n)
        const arma::mat & Theta,            // (p,d)
        const arma::mat & B,                // (p,q)
        const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,q)

    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &Theta, &B](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M * B.t();
        arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
        arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2) ;
        double objective = accu(diagmat(w) * (A - Y % Z)) + 0.5 * accu(diagmat(w) * (M % M + S2 - log(S2) - 1.));

        metadata.map<M_ID>(grad) = diagmat(w) * ((A - Y) * B + M);
        metadata.map<S_ID>(grad) = diagmat(w) * (S - 1. / S + A * (B % B) % S);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M * B.t();
    arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
    arma::mat loglik = arma::sum(Y % Z - A, 1) - 0.5 * sum(M % M + S2 - log(S2) - 1., 1) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
}
