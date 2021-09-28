#include <RcppArmadillo.h>

#include "nlopt_wrapper.h"
#include "packing.h"

// ---------------------------------------------------------------------------------------
// Step a, returns new Omega (p,p)
// could be R only

// [[Rcpp::export]]
arma::mat cpp_optimize_zi_Omega(
    const arma::mat & M,     // (n,p)
    const arma::mat & X,     // (n,d)
    const arma::mat & Theta, // (d,p)
    const arma::mat & S     // s_{i,j}^2 (n,p)
) {
    const arma::uword n = M.n_rows;
    arma::mat M_X_Theta = M - X * Theta;
    return (1. / double(n)) * inv_sympd(M_X_Theta.t() * M_X_Theta + diagmat(sum(S % S, 0)));
}

// ---------------------------------------------------------------------------------------
// Step b, returns new theta (d,p)
// could be R only

// [[Rcpp::export]]
arma::mat cpp_optimize_zi_Theta(
    const arma::mat & M, // (n,p)
    const arma::mat & X  // (n,d)
) {
    // Armadillo requires using solve(A,B) for Ax=B <=> x = A^-1 B
    // So (X^T X)^-1 X^T M becomes solve(A,B) with A = X^T X and B = X^T M.
    // X^T X is sympd, provide this indications to solve()
    return solve(X.t() * X, X.t() * M, arma::solve_opts::likely_sympd);
}

// ---------------------------------------------------------------------------------------
// Step c, optimizes theta0

// [[Rcpp::export]]
Rcpp::List cpp_optimize_zi_Theta0(
    const arma::mat & init_Theta0,   // (d,p)
    const arma::mat & X,             // covariates (n,d)
    const arma::mat & Pi,            // (n,p)
    const Rcpp::List & configuration // List of config values ; xtol_abs is Theta0 only (double or mat)
) {
    const auto metadata = tuple_metadata(init_Theta0);
    enum { THETA0_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA0_ID>(parameters.data()) = init_Theta0;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto packed = std::vector<double>(metadata.packed_size);
            metadata.map<THETA0_ID>(packed.data()) = Rcpp::as<arma::mat>(value);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const arma::mat Xt_Pi = X.t() * Pi;

    // Optimize
    auto objective_and_grad = [&metadata, &X, &Pi, &Xt_Pi](const double * params, double * grad) -> double {
        const arma::mat Theta0 = metadata.map<THETA0_ID>(params);

        arma::mat e_X_Theta0 = exp(X * Theta0);
        double objective = -trace(Xt_Pi.t() * Theta0) + accu(log(1. + e_X_Theta0));
        metadata.map<THETA0_ID>(grad) = -Xt_Pi + X.t() * pow(1. + e_X_Theta0, -1);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat Theta0 = metadata.copy<THETA0_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("Theta0") = Theta0);
}

// ---------------------------------------------------------------------------------------
// Step d, returns new Pi (n,p)

// [[Rcpp::export]]
arma::mat cpp_optimize_zi_Pi(
    const arma::mat & Y,     // responses (n,p)
    const arma::mat & X,     // covariates (n,d)
    const arma::mat & O,     // offsets (n,p)
    const arma::mat & M,     // (n,p)
    const arma::mat & S,    // (n,p)
    const arma::mat & Theta0 // (d,p)
) {
    arma::mat A = exp(O + M + 0.5 * S % S);
    arma::mat Pi = pow(1. + exp(- (A + X * Theta0)), -1);
    // Zero Pi_{i,j} if Y_{i,j} > 0
    // multiplication with f(sign(Y)) could work to zero stuff as there should not be any +inf
    // using a loop as it is more explicit and should have ok performance in C++
    arma::uword n = Y.n_rows;
    arma::uword p = Y.n_cols;
    for(arma::uword i = 0; i < n; i += 1) {
        for(arma::uword j = 0; j < p; j += 1) {
            // Add fuzzy comparison ?
            if(Y(i, j) > 0.) {
                Pi(i, j) = 0.;
            }
        }
    }
    return Pi;
}

// ---------------------------------------------------------------------------------------
// Step e, optimizes M

// [[Rcpp::export]]
Rcpp::List cpp_optimize_zi_M(
    const arma::mat & init_M,        // (n,p)
    const arma::mat & Y,             // responses (n,p)
    const arma::mat & X,             // covariates (n,d)
    const arma::mat & O,             // offsets (n, p)
    const arma::mat & Pi,            // (n,p)
    const arma::mat & S,             // (n,p)
    const arma::mat & Theta,         // (d,p)
    const arma::mat & Omega,         // (p,p)
    const Rcpp::List & configuration // List of config values ; xtol_abs is M only (double or mat)
) {
    const auto metadata = tuple_metadata(init_M);
    enum { M_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto packed = std::vector<double>(metadata.packed_size);
            metadata.map<M_ID>(packed.data()) = Rcpp::as<arma::mat>(value);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const arma::mat X_Theta = X * Theta; // (n,p)
    const arma::mat O_S2 = O + 0.5 * S % S; // (n,p)

    // Optimize
    auto objective_and_grad =
        [&metadata, &Y, &X, &O_S2, &Pi, &X_Theta, &Omega](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);

        arma::mat A = exp(O_S2 + M);                       // (n,p)
        arma::mat M_X_Theta_Omega = (M - X_Theta) * Omega; // (n,p)

        double objective = trace((Pi - 1.).t() * (Y % M - A)) + 0.5 * trace(M_X_Theta_Omega * (M - X_Theta).t());
        metadata.map<M_ID>(grad) = (Pi - 1.) % (Y - A) + M_X_Theta_Omega;
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M = metadata.copy<M_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M);
}

// ---------------------------------------------------------------------------------------
// Step f, optimizes S

// [[Rcpp::export]]
Rcpp::List cpp_optimize_zi_S(
    const arma::mat & init_S,        // (n,p)
    const arma::mat & O,             // offsets (n, p)
    const arma::mat & M,             // (n,p)
    const arma::mat & Pi,            // (n,p)
    const arma::mat & Theta,         // (d,p)
    const arma::mat & Omega,         // (p,p)
    const Rcpp::List & configuration // List of config values ; xtol_abs is S2 only (double or mat)
) {
    const auto metadata = tuple_metadata(init_S);
    enum { S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto packed = std::vector<double>(metadata.packed_size);
            metadata.map<S_ID>(packed.data()) = Rcpp::as<arma::mat>(value);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const arma::mat O_M = O + M;
    const arma::vec diag_Omega = diagvec(Omega);

    // Optimize
    auto objective_and_grad = [&metadata, &O_M, &Pi, &diag_Omega](const double * params, double * grad) -> double {
        const arma::mat S = metadata.map<S_ID>(params);

        arma::uword n = S.n_rows;
        arma::mat A = exp(O_M + 0.5 * S % S); // (n,p)

        // trace(1^T log(S)) == accu(log(S)).
        // S_bar = diag(sum(S, 0)). trace(Omega * S_bar) = dot(diagvec(Omega), sum(S2, 0))
        double objective = trace((Pi - 1.).t() * A) + 0.5 * dot(diag_Omega, sum(S % S, 0)) - 0.5 * accu(log(S % S));
        // S2^\emptyset interpreted as pow(S2, -1.) as that makes the most sense (gradient component for log(S2))
        // 1_n Diag(Omega)^T is n rows of diag(omega) values
        metadata.map<S_ID>(grad) = - pow(S, -1.) + (1. - Pi) % S % A + S * diagmat(diag_Omega);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat S = metadata.copy<S_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("S") = S);
}
