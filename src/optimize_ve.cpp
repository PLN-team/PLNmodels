#include <RcppArmadillo.h>

#include "nlopt_wrapper.h"
#include "packer.h"

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
    return -logfact(std::move(y)) + 0.5 * (1. + (1. - double(p)) * std::log(2. * M_PI));
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,p)

    const auto packer = make_packer(init_M, init_S);
    enum { M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &Theta, &Omega](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat nSigma = M.t() * diagmat(w) * M + diagmat(sum(S2.each_col() % w, 0));
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * trace(Omega * nSigma);

        packer.pack<M_ID>(grad_storage, diagmat(w) * (M * Omega + A - Y));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,p)

    const auto packer = make_packer(init_M, init_S);
    enum { M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &Theta, &Omega](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S);
        arma::vec omega2 = arma::diagvec(Omega);
        double objective =
            accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * as_scalar(w.t() * (pow(M, 2) + S2) * omega2);

        packer.pack<M_ID>(grad_storage, diagmat(w) * (M * Omega + A - Y));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S.each_row() % omega2 + S % A - pow(S, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::vec>(init_parameters["S"]); // (n)

    const auto packer = make_packer(init_M, init_S);
    enum { M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &Theta, &Omega](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::vec S = packer.unpack<S_ID>(parameters);

        arma::vec S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z.each_col() + 0.5 * S2);
        const arma::uword p = Y.n_cols;
        double n_sigma2 = dot(w, sum(pow(M, 2), 1) + double(p) * S);
        double omega2 = Omega(0, 0);
        double objective = accu(w.t() * (A - Y % Z)) - 0.5 * double(p) * dot(w, log(S2)) + 0.5 * n_sigma2 * omega2;

        packer.pack<M_ID>(grad_storage, diagmat(w) * (M * omega2 + A - Y));
        packer.pack<S_ID>(grad_storage, w % (S % sum(A, 1) - double(p) * pow(S, -1) - double(p) * S * omega2));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters); // vec(n) -> mat(n, 1)
    arma::vec S2 = S % S;
    double omega2 = Omega(0, 0);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z.each_col() + 0.5 * S2);
    const arma::uword p = Y.n_cols;
    arma::mat loglik = sum(Y % Z - A - 0.5 * pow(M, 2) * omega2, 1) - 0.5 * double(p) * omega2 * S2 +
                       0.5 * double(p) * log(S2 * omega2) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status") = (int)result.status,
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S,
        Rcpp::Named("loglik") = loglik);
}
