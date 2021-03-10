#include <RcppArmadillo.h>

#include "nlopt_wrapper.h"
#include "packing.h"

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
//    return -logfact(std::move(y)) + 0.5 * (1. + (1. - double(p)) * std::log(2. * M_PI));
    return -logfact(std::move(y)) + 0.5 * double(p) ;
}
// ---------------------------------------------------------------------------------------
// VE full

// [[Rcpp::export]]
Rcpp::List cpp_optimize_vestep_full(
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
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,p)

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
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2) ;
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * trace(Omega * nSigma);

        metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A + 0.5 * log(S2) - 0.5 * ((M * Omega) % M + S * diagmat(Omega)), 1) +
                       0.5 * real(log_det(Omega)) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
}

// ---------------------------------------------------------------------------------------
// VE diagonal

// [[Rcpp::export]]
Rcpp::List cpp_optimize_vestep_diagonal(
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
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,p)

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
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::vec omega2 = arma::diagvec(Omega);
        double objective =
            accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * as_scalar(w.t() * (pow(M, 2) + S2) * omega2) ;

        metadata.map<M_ID>(grad) = diagmat(w) * (M * arma::diagmat(omega2) + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % omega2.t() + S % A - pow(S, -1));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::vec omega2 = Omega.diag();
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik =
        sum(Y % Z - A + 0.5 * log(S2), 1) - 0.5 * (pow(M, 2) + S2) * omega2 + 0.5 * sum(log(omega2)) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
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
    const auto init_S = Rcpp::as<arma::vec>(init_parameters["S"]); // (n)

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

        arma::vec S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2 * arma::ones(p).t());
        double n_sigma2 = dot(w, sum(pow(M, 2), 1) + double(p) * S2);
        double omega2 = Omega(0, 0);
        double objective = accu(w.t() * (A - Y % Z)) - 0.5 * double(p) * dot(w, log(S2)) + 0.5 * n_sigma2 * omega2;

        metadata.map<M_ID>(grad) = diagmat(w) * (M * omega2 + A - Y);
        metadata.map<S_ID>(grad) = w % (double(p) * S * omega2 + S % sum(A, 1) - double(p) * pow(S, -1));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data()); // vec(n) -> mat(n, 1)
    arma::vec S2 = S % S;
    double omega2 = Omega(0, 0);
    // Element-wise log-likelihood
    const arma::uword p = Y.n_cols;
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2 * arma::ones(p).t());
    arma::mat loglik = sum(Y % Z - A - 0.5 * pow(M, 2) * omega2, 1) - 0.5 * double(p) * omega2 * S2 +
                       0.5 * double(p) * log(S2 * omega2) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
}
