#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "lambertW.h"

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
    const arma::mat & S       // (n,p)
) {
    const arma::uword p = Y.n_cols;

    const arma::mat S2 = S % S ;
    const arma::mat A = exp(O + M + .5 * S2) ;
    const arma::mat M_mu = M - X * B ;
    const arma::mat mu0  = logit(Pi) ;
    return (
        0.5 * real(log_det(Omega)) + 0.5 * double(p)
        + sum(
            (1 - R) % ( Y % (O + M) - A - logfact_mat(Y) )
            + R % mu0 - log( 1 + exp(mu0) )
            + 0.5 * log(S2) - 0.5 * ((M_mu * Omega) % M_mu + S2 * diagmat(Omega))
            - R % trunc_log(R) - (1 - R) % trunc_log(1-R), 1)
    ) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_full(
    const arma::mat & M, // (n,p)
    const arma::mat & X, // (n,d)
    const arma::mat & B, // (d,p)
    const arma::mat & S  // (n,p)
) {
    const arma::uword n = M.n_rows;
    arma::mat M_mu = M - X * B;
    return (double(n) * inv_sympd(M_mu.t() * M_mu + diagmat(sum(S % S, 0))));
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_spherical(
    const arma::mat & M, // (n,p)
    const arma::mat & X, // (n,d)
    const arma::mat & B, // (d,p)
    const arma::mat & S  //  (n,p)
) {
    const arma::uword n = M.n_rows;
    const arma::uword p = M.n_cols;
    double sigma2 = accu( pow(M - X * B, 2) + S % S ) / double(n * p) ;
    return arma::diagmat(arma::ones(p)/sigma2) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_Omega_diagonal(
    const arma::mat & M, // (n,p)
    const arma::mat & X, // (n,d)
    const arma::mat & B, // (d,p)
    const arma::mat & S  // (n,p)
) {
    const arma::uword n = M.n_rows;
    return arma::diagmat(double(n) / sum( pow(M - X * B, 2) + S % S, 0)) ;
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

    // Optimize
    auto objective_and_grad = [&metadata, &X0, &R, &Xt_R](const double * params, double * grad) -> double {
        const arma::mat B0 = metadata.map<B0_ID>(params);

        arma::mat e_mu0 = exp(X0 * B0);
        double objective = -trace(Xt_R.t() * B0) + accu(log(1. + e_mu0));
        metadata.map<B0_ID>(grad) = -Xt_R + X0.t() * (e_mu0 % pow(1. + e_mu0, -1)) ;
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
arma::mat optim_zipln_R_var(
    const arma::mat & Y, // responses (n,p)
    const arma::mat & X, // covariates (n,d)
    const arma::mat & O, // offsets (n,p)
    const arma::mat & M, // (n,p)
    const arma::mat & S, // (n,p)
    const arma::mat & Pi, // (d,p)
    const arma::mat & B // covariates (n,d)
) {
    arma::mat A = exp(O + M + 0.5 * S % S);
    arma::mat R = pow(1. + exp(- (A + logit(Pi))), -1);
    // Zero R_{i,j} if Y_{i,j} > 0
    // multiplication with f(sign(Y)) could work to zero stuff as there should not be any +inf
    // using a loop as it is more explicit and should have ok performance in C++
    arma::uword n = Y.n_rows;
    arma::uword p = Y.n_cols;
    for(arma::uword i = 0; i < n; i += 1) {
        for(arma::uword j = 0; j < p; j += 1) {
            // Add fuzzy comparison ?
            if(Y(i, j) > 0.) {
                R(i, j) = 0.;
            }
        }
    }
    return R;
}

double phi (double mu, double sigma2) {
  double W = lambertW0_CS(sigma2 * exp(mu)) ;
  return(exp(-(pow(W, 2) + 2 * W) / (2 * sigma2)) / sqrt(1 + W)) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_R_exact (
    const arma::mat & Y, // covariates (n,d)
    const arma::mat & X, // covariates (n,d)
    const arma::mat & O, // offsets (n,p)
    const arma::mat & M, // (n,p)
    const arma::mat & S, // (n,p)
    const arma::mat & Pi, // (n,p)
    const arma::mat & B // covariates (n,d)
) {

  arma::mat XB = X * B;
  arma::mat M_mu = M - XB;
  arma::uword n = M.n_rows;
  arma::uword p = M.n_cols;
  arma::vec diag_Sigma = arma::diagvec((1./n) * (M_mu.t() * M_mu + diagmat(sum(S % S, 0)))) ;
  arma::mat R = arma::zeros(n,p);
  for(arma::uword i = 0; i < n; i += 1) {
    for(arma::uword j = 0; j < p; j += 1) {
      if(Y(i, j) < 0.5) {
        double Phi = phi(O(i,j) + XB(i,j), diag_Sigma(j)) ;
        R(i,j) = Pi(i,j) / (Phi * (1 - Pi(i,j)) + Pi(i,j)) ;
      }
    }
  }
  return R;
}

// [[Rcpp::export]]
Rcpp::List optim_zipln_M(
    const arma::mat & init_M, // (n,p)
    const arma::mat & Y,      // responses (n,p)
    const arma::mat & X,      // covariates (n,d)
    const arma::mat & O,      // offsets (n, p)
    const arma::mat & R,      // (n,p)
    const arma::mat & S,      // (n,p)
    const arma::mat & B,      // (d,p)
    const arma::mat & Omega,  // (p,p)
    const Rcpp::List & configuration // List of config values ; xtol_abs is M only (double or mat)
) {
    const auto metadata = tuple_metadata(init_M);
    enum { M_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat X_B = X * B; // (n,p)
    const arma::mat O_S2 = O + 0.5 * S % S; // (n,p)

    // Optimize
    auto objective_and_grad =
        [&metadata, &Y, &X, &O_S2, &R, &X_B, &Omega](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);

        arma::mat A = exp(O_S2 + M);              // (n,p)
        arma::mat M_mu_Omega = (M - X_B) * Omega; // (n,p)

        double objective = - trace((1. - R).t() * (Y % M - A)) + 0.5 * trace(M_mu_Omega * (M - X_B).t());
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

// [[Rcpp::export]]
Rcpp::List optim_zipln_S(
    const arma::mat & init_S,    // (n,p)
    const arma::mat & O,         // offsets (n, p)
    const arma::mat & M,         // (n,p)
    const arma::mat & R,         // (n,p)
    const arma::mat & B,         // (d,p)
    const arma::vec & diag_Omega,// (p,1)
    const Rcpp::List & configuration // List of config values ; xtol_abs is S2 only (double or mat)
) {
    const auto metadata = tuple_metadata(init_S);
    enum { S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat O_M = O + M;

    // Optimize
    auto objective_and_grad = [&metadata, &O_M, &R, &diag_Omega](const double * params, double * grad) -> double {
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat A = exp(O_M + 0.5 * S % S); // (n,p)

        // trace(1^T log(S)) == accu(log(S)).
        // S_bar = diag(sum(S, 0)). trace(Omega * S_bar) = dot(diagvec(Omega), sum(S2, 0))
        double objective = trace((1. - R).t() * A) + 0.5 * dot(diag_Omega, sum(S % S, 0)) - 0.5 * accu(log(S % S));
        // S2^\emptyset interpreted as pow(S2, -1.) as that makes the most sense (gradient component for log(S2))
        // 1_n Diag(Omega)^T is n rows of diag(omega) values
        metadata.map<S_ID>(grad) = S.each_row() % diag_Omega.t() + (1. - R) % S % A - pow(S, -1.) ;
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat S = metadata.copy<S_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("S") = S);
}
