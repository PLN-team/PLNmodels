#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]


#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Spherical covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_spherical(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
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
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &w_bar](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * (double(p) * w_bar) * log(sigma2) ;

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * (M / sigma2 + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S / sigma2 + S % A - pow(S, -1));

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    // Variance parameters
    const arma::uword p = Y.n_cols;
    const double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
    arma::sp_mat Sigma(p,p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p,p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2 ) / sigma2 + 0.5 * log(S2 / sigma2), 1) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", Theta),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// VE spherical

// [[Rcpp::export]]
Rcpp::List cpp_optimize_vestep_spherical(
    const Rcpp::List & init_parameters, // List(M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const arma::mat & Theta,            // (p,d)
    const arma::mat & Omega,            // (p,p)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n)

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
    auto objective_and_grad =
        [&metadata, &O, &X, &Y, &w, &Theta, &Omega](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double n_sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) ;
        double omega2 = Omega(0, 0);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * n_sigma2 * omega2;

        metadata.map<M_ID>(grad) = diagmat(w) * (M / omega2 + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S / omega2 + S % A - pow(S, -1));

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    double omega2 = Omega(0, 0);
    // Element-wise log-likelihood
    const arma::uword p = Y.n_cols;
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2 ) * omega2 + 0.5 * log(S2 * omega2), 1) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
}
