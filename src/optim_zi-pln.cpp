#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// All optimizers use the reparameterization ψ = log(S²) so the variance S² is
// always positive. The R/C++ interface passes S2 (variance) rather than S (std dev).

// [[Rcpp::export]]
arma::vec zipln_vloglik(
    const arma::mat & Y,      // responses (n,p)
    const arma::mat & X,      // covariates (n,d)
    const arma::mat & O,      // offsets (n,p)
    const arma::mat & Pi,     // (d,p)
    const arma::mat & Omega,  // (p,p)
    const arma::mat & B,      // (d,p)
    const arma::mat & R,      // (n,p)
    const arma::mat & M,      // (n,p)
    const arma::mat & S2      // (n,p) variational variance
) {
    const arma::uword p = Y.n_cols;

    const arma::mat A    = trunc_exp(O + M + .5 * S2) ;
    const arma::mat M_mu = M - X * B ;
    const arma::mat mu0  = logit(Pi) ;
    return (
        0.5 * real(log_det(Omega)) + 0.5 * double(p)
        + sum(
            (1 - R) % ( Y % (O + M) - A - logfact_mat(Y) )
            + R % mu0 - trunc_log( 1 + exp(mu0) )
            + 0.5 * trunc_log(S2) - 0.5 * ((M_mu * Omega) % M_mu + S2.each_row() % diagvec(Omega).t())
            - R % trunc_log(R) - (1 - R) % trunc_log(1-R), 1)
    ) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_full(
    const arma::mat & M,  // (n,p)
    const arma::mat & X,  // (n,d)
    const arma::mat & B,  // (d,p)
    const arma::mat & S2  // (n,p) variational variance
) {
    const arma::uword n = M.n_rows;
    arma::mat M_mu = M - X * B;
    return (double(n) * inv_sympd(M_mu.t() * M_mu + diagmat(sum(S2, 0))));
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_spherical(
    const arma::mat & M,  // (n,p)
    const arma::mat & X,  // (n,d)
    const arma::mat & B,  // (d,p)
    const arma::mat & S2  // (n,p) variational variance
) {
    const arma::uword n = M.n_rows;
    const arma::uword p = M.n_cols;
    double sigma2 = accu( arma::square(M - X * B) + S2 ) / double(n * p) ;
    return arma::diagmat(arma::ones(p)/sigma2) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_diagonal(
    const arma::mat & M,  // (n,p)
    const arma::mat & X,  // (n,d)
    const arma::mat & B,  // (d,p)
    const arma::mat & S2  // (n,p) variational variance
) {
    const arma::uword n = M.n_rows;
    return arma::diagmat(double(n) / sum( arma::square(M - X * B) + S2, 0)) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_B_dense(
    const arma::mat & M, // (n,p)
    const arma::mat & X  // (n,d)
) {
    // X^T X is sympd, provide this indications to solve()
    return solve(X.t() * X, X.t() * M, arma::solve_opts::likely_sympd);
}

// [[Rcpp::export]]
Rcpp::List optim_zipln_zipar_covar(
    const arma::mat & R,        // (n,p)
    const arma::mat & init_B0,  // (d0,p)
    const arma::mat & X0,       // covariates (n,d0)
    const Rcpp::List & configuration // List of config values ; xtol_abs is B0 only (double or mat)
) {
    const auto metadata = tuple_metadata(init_B0);
    enum { B0_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B0_ID>(parameters.data()) = init_B0;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());

    const arma::mat Xt_R = X0.t() * R;
    const arma::mat X0t = X0.t();

    // Optimize
    auto objective_and_grad = [&metadata, &X0, &X0t, &Xt_R](const double * params, double * grad) -> double {
        const arma::mat B0 = metadata.map<B0_ID>(params);

        arma::mat e_mu0 = exp(X0 * B0);
        double objective = -accu(Xt_R % B0) + accu(log(1. + e_mu0));
        metadata.map<B0_ID>(grad) = -Xt_R + X0t * (e_mu0 / (1. + e_mu0)) ;
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat B0 = metadata.copy<B0_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("B0") = B0,
        Rcpp::Named("Pi") = logistic(X0 * B0));
}

// [[Rcpp::export]]
Rcpp::List optim_zipln_M(
    const arma::mat & init_M, // (n,p)
    const arma::mat & Y,      // responses (n,p)
    const arma::mat & X,      // covariates (n,d)
    const arma::mat & O,      // offsets (n, p)
    const arma::mat & R,      // (n,p)
    const arma::mat & S2,     // (n,p) variational variance (fixed)
    const arma::mat & B,      // (d,p)
    const arma::mat & Omega,  // (p,p)
    const Rcpp::List & configuration // List of config values ; xtol_abs is M only (double or mat)
) {
    const auto metadata = tuple_metadata(init_M);
    enum { M_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat X_B  = X * B;           // (n,p)
    const arma::mat O_S2 = O + 0.5 * S2;   // (n,p)

    // Optimize
    auto objective_and_grad =
        [&metadata, &Y, &X, &O_S2, &R, &X_B, &Omega](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);

        arma::mat A = exp(O_S2 + M);              // (n,p)
        arma::mat M_mu_Omega = (M - X_B) * Omega; // (n,p)

        double objective = - accu((1. - R) % (Y % M - A)) + 0.5 * accu(M_mu_Omega % (M - X_B));
        metadata.map<M_ID>(grad) = M_mu_Omega + (1. - R) % (A - Y);
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
// Optimize ψ = log(S²) only, M fixed — nlopt/CCSAQ
// Interface: takes S2 (variance), returns S2.

// [[Rcpp::export]]
Rcpp::List optim_zipln_psi(
    const arma::mat & init_S2,   // (n,p) variational variance (initialization)
    const arma::mat & O,         // offsets (n, p)
    const arma::mat & M,         // (n,p) fixed
    const arma::mat & R,         // (n,p)
    const arma::mat & B,         // (d,p)
    const arma::vec & diag_Omega,// (p)
    const Rcpp::List & configuration
) {
    const arma::mat psi_init = arma::log(init_S2);
    const auto metadata = tuple_metadata(psi_init);
    enum { PSI_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<PSI_ID>(parameters.data()) = psi_init;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat O_M = O + M;

    auto objective_and_grad = [&metadata, &O_M, &R, &diag_Omega](const double * params, double * grad) -> double {
        const arma::mat psi = metadata.map<PSI_ID>(params);
        const arma::mat S2  = arma::exp(psi);
        const arma::mat A   = exp(O_M + 0.5 * S2);

        // f = accu((1-R)%A) + 0.5*dot(diag_Omega, sum(S2,0)) - 0.5*accu(psi)
        double objective = accu((1. - R) % A) + 0.5 * dot(diag_Omega, sum(S2, 0)) - 0.5 * accu(psi);

        // grad_ψ_ij = 0.5 * S2_ij * (diag_Omega_j + (1-R_ij)*A_ij) - 0.5
        metadata.map<PSI_ID>(grad) = 0.5 * (S2.each_row() % diag_Omega.t() + (1. - R) % S2 % A - 1.);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat psi = metadata.copy<PSI_ID>(parameters.data());
    arma::mat S2  = arma::exp(psi);
    return Rcpp::List::create(
        Rcpp::Named("status")     = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("S2") = S2);
}

// VE step for (M, ψ=log(S²), R) — ported to the generic, trait-based
// ve_step_zipln_{newton,nlopt}_impl<Traits> templates (builtin_optim_zipln.h /
// nlopt_optim_zipln.h), exported per covariance structure in
// wrappers_builtin_optim_zipln.cpp / wrappers_nlopt_optim_zipln.cpp.
