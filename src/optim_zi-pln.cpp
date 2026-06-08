#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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
    const arma::mat A = trunc_exp(O + M + .5 * S2) ;
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
    arma::mat R = 1. / (1. + exp(-(A + logit(Pi))));
    // Zero R_{i,j} if Y_{i,j} > 0
    // multiplication with f(sign(Y)) could work to zero stuff as there should not be any +inf
    // using a loop as it is more explicit and should have ok performance in C++
    R.elem(arma::find(Y > 0.)).zeros();
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
  const int n = (int)M.n_rows;
  const int p = (int)M.n_cols;
  arma::vec diag_Sigma = (sum(M_mu % M_mu, 0) + sum(S % S, 0)).t() / double(n);
  arma::mat R = arma::zeros(n, p);
  // lambertW0_CS is pure (no global state) — safe to parallelize
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < p; j++) {
      if(Y(i, j) < 0.5) {
        double Phi = phi(O(i,j) + XB(i,j), diag_Sigma(j));
        R(i,j) = Pi(i,j) / (Phi * (1 - Pi(i,j)) + Pi(i,j));
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

        const arma::mat S2 = S % S;
        arma::mat A = exp(O_M + 0.5 * S2); // (n,p)

        // trace(1^T log(S)) == accu(log(S)).
        // S_bar = diag(sum(S, 0)). trace(Omega * S_bar) = dot(diagvec(Omega), sum(S2, 0))
        double objective = accu((1. - R) % A) + 0.5 * dot(diag_Omega, sum(S2, 0)) - 0.5 * accu(log(S2));
        // S2^\emptyset interpreted as pow(S2, -1.) as that makes the most sense (gradient component for log(S2))
        // 1_n Diag(Omega)^T is n rows of diag(omega) values
        metadata.map<S_ID>(grad) = S.each_row() % diag_Omega.t() + (1. - R) % S % A - 1. / S ;
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat S = metadata.copy<S_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status") = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("S") = S);
}

// [[Rcpp::export]]
Rcpp::List optim_zipln_M_S(
    const arma::mat & init_M,    // (n,p)
    const arma::mat & init_S,    // (n,p)
    const arma::mat & Y,         // responses (n,p)
    const arma::mat & X,         // covariates (n,d)
    const arma::mat & O,         // offsets (n,p)
    const arma::mat & R,         // (n,p)
    const arma::mat & B,         // (d,p)
    const arma::mat & Omega,     // (p,p)
    const Rcpp::List & configuration
) {
    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat X_B = X * B;
    const arma::vec diag_Omega = diagvec(Omega);

    auto objective_and_grad = [&metadata, &Y, &O, &R, &X_B, &Omega, &diag_Omega](
            const double * params, double * grad) -> double {
        const arma::mat M  = metadata.map<M_ID>(params);
        const arma::mat S  = metadata.map<S_ID>(params);
        const arma::mat S2 = S % S;
        const arma::mat A  = exp(O + M + 0.5 * S2);
        const arma::mat M_mu        = M - X_B;
        const arma::mat M_mu_Omega  = M_mu * Omega;

        double objective = - accu((1. - R) % (Y % M - A))
                         + 0.5 * accu(M_mu_Omega % M_mu)
                         + 0.5 * dot(diag_Omega, sum(S2, 0))
                         - 0.5 * accu(log(S2));

        metadata.map<M_ID>(grad) = M_mu_Omega + (1. - R) % (A - Y);
        metadata.map<S_ID>(grad) = S.each_row() % diag_Omega.t() + (1. - R) % S % A - 1. / S;

        return objective;
    };

    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status")     = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = S
    );
}

// ---------------------------------------------------------------------------------------
// Joint optimization of (M, logS) — reparametrize S = exp(logS) to make the problem
// unconstrained, enabling quasi-Newton methods (LBFGS) that diverge in the S domain.
//
// Gradient w.r.t. logS_ij = S_ij * grad_S_ij
//   = S²_ij * (diag_Omega_j + (1-R_ij)*A_ij) - 1
// Hessian diagonal (exact, S decoupled):
//   = 2*S²_ij * (diag_Omega_j + (1-R_ij)*A_ij*(1 + S²_ij))  > 0
// This makes the problem globally unconstrained with positive-definite diagonal Hessian in logS.

// [[Rcpp::export]]
Rcpp::List optim_zipln_M_logS(
    const arma::mat & init_M,    // (n,p)
    const arma::mat & init_S,    // (n,p)  — converted to logS internally
    const arma::mat & Y,         // (n,p)
    const arma::mat & X,         // (n,d)
    const arma::mat & O,         // (n,p)
    const arma::mat & R,         // (n,p)
    const arma::mat & B,         // (d,p)
    const arma::mat & Omega,     // (p,p)
    const Rcpp::List & configuration
) {
    const arma::mat logS_init = arma::log(init_S);

    const auto metadata = tuple_metadata(init_M, logS_init);
    enum { M_ID, logS_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>  (parameters.data()) = init_M;
    metadata.map<logS_ID>(parameters.data()) = logS_init;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat X_B = X * B;
    const arma::vec diag_Omega = diagvec(Omega);

    auto objective_and_grad = [&metadata, &Y, &O, &R, &X_B, &Omega, &diag_Omega](
            const double * params, double * grad) -> double {
        const arma::mat M    = metadata.map<M_ID>  (params);
        const arma::mat logS = metadata.map<logS_ID>(params);
        const arma::mat S2   = arma::exp(2. * logS);
        const arma::mat A    = arma::exp(O + M + 0.5 * S2);
        const arma::mat M_mu       = M - X_B;
        const arma::mat M_mu_Omega = M_mu * Omega;

        double objective = - accu((1. - R) % (Y % M - A))
                         + 0.5 * accu(M_mu_Omega % M_mu)
                         + 0.5 * dot(diag_Omega, sum(S2, 0))
                         - accu(logS);   // = -0.5*accu(log(S²))

        metadata.map<M_ID>  (grad) = M_mu_Omega + (1. - R) % (A - Y);
        // grad_logS_ij = S²_ij * (diag_Omega_j + (1-R_ij)*A_ij) - 1
        metadata.map<logS_ID>(grad) = S2 % ((1. - R) % A + arma::ones(A.n_rows, 1) * diag_Omega.t()) - 1.;

        return objective;
    };

    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M    = metadata.copy<M_ID>  (parameters.data());
    arma::mat logS = metadata.copy<logS_ID>(parameters.data());
    return Rcpp::List::create(
        Rcpp::Named("status")     = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = arma::exp(logS)
    );
}

// ---------------------------------------------------------------------------------------
// Custom coordinate-Newton optimizer for (M, logS) — no NLopt dependency.
//
// Per iteration:
//   1. Diagonal Newton step for M with Armijo backtracking (one M*Omega product)
//   2. Exact closed-form update for logS: logS* = -0.5*log(diag_Omega + (1-R)*A)
//      This is the exact minimiser of f w.r.t. S for fixed M and A (before S feeds back
//      into A). Stable by construction; no step-size issue.
//
// Cost per iteration ≈ 2 × O(n*p²) gradient evaluations — same arithmetic as CCSAQ
// but with much better step quality: typically 5-15 iterations vs 50-100 for CCSAQ.

// [[Rcpp::export]]
Rcpp::List optim_zipln_M_S_newton(
    const arma::mat & init_M,    // (n,p)
    const arma::mat & init_S,    // (n,p)
    const arma::mat & Y,         // (n,p)
    const arma::mat & X,         // (n,d)
    const arma::mat & O,         // (n,p)
    const arma::mat & R,         // (n,p)
    const arma::mat & B,         // (d,p)
    const arma::mat & Omega,     // (p,p)
    const int    maxiter,        // max coordinate steps
    const double ftol_rel        // relative objective convergence threshold
) {
    const arma::mat X_B        = X * B;
    const arma::vec diag_Omega = diagvec(Omega);
    const arma::mat ones_col   = arma::ones(init_M.n_rows, 1);

    arma::mat M    = init_M;
    arma::mat logS = arma::log(init_S);

    // Helper: full objective value (for convergence check and Armijo)
    auto objective = [&](const arma::mat & M_, const arma::mat & logS_) -> double {
        const arma::mat S2_   = arma::exp(2. * logS_);
        const arma::mat A_    = arma::exp(O + M_ + 0.5 * S2_);
        const arma::mat M_mu_ = M_ - X_B;
        return - accu((1. - R) % (Y % M_ - A_))
               + 0.5 * accu((M_mu_ * Omega) % M_mu_)
               + 0.5 * dot(diag_Omega, sum(S2_, 0))
               - accu(logS_);
    };

    double obj_prev = objective(M, logS);
    int iter = 0;

    for (iter = 0; iter < maxiter; iter++) {
        // ----------------------------------------------------------------
        // Step 1 — diagonal Newton on M with Armijo backtracking
        // ----------------------------------------------------------------
        {
            const arma::mat S2         = arma::exp(2. * logS);
            const arma::mat A          = arma::exp(O + M + 0.5 * S2);
            const arma::mat M_mu_Omega = (M - X_B) * Omega;
            const arma::mat grad_M     = M_mu_Omega + (1. - R) % (A - Y);
            const arma::mat hess_M     = (1. - R) % A + ones_col * diag_Omega.t();
            const arma::mat step_M     = grad_M / hess_M;   // Newton direction (descent)

            // Armijo backtracking: f(M - alpha*step) <= f(M) - c1*alpha*||step||²_H
            double alpha   = 1.0;
            const double c1 = 1e-4;
            const double slope = - accu(grad_M % step_M);   // ≤ 0
            for (int ls = 0; ls < 20; ++ls) {
                if (objective(M - alpha * step_M, logS) <= obj_prev + c1 * alpha * slope)
                    break;
                alpha *= 0.5;
            }
            M -= alpha * step_M;
        }

        // ----------------------------------------------------------------
        // Step 2 — exact closed-form update for logS (fixed-point step)
        //   logS* = -0.5 * log( diag_Omega + (1-R)*A )
        // Per-element upper bound prevents exp(Z + 0.5*S²) from overflowing:
        //   keep Z + 0.5*S² ≤ 700 => logS ≤ 0.5*log(max(1, 700-Z))
        // ----------------------------------------------------------------
        {
            const arma::mat S2       = arma::exp(2. * logS);
            const arma::mat Z_cur    = O + X_B + M;
            const arma::mat A        = arma::exp(Z_cur + 0.5 * S2);
            const arma::mat base     = (1. - R) % A + ones_col * diag_Omega.t();
            const arma::mat logS_cand = -0.5 * arma::log(base);
            const arma::mat logS_ub  = 0.5 * arma::log(
                arma::clamp(700. - Z_cur, 1., arma::datum::inf));
            logS = arma::clamp(arma::min(logS_cand, logS_ub), -20., arma::datum::inf);
        }

        // ----------------------------------------------------------------
        // Convergence check
        // ----------------------------------------------------------------
        const double obj = objective(M, logS);
        if (iter > 0 && std::abs(obj - obj_prev) < ftol_rel * (1. + std::abs(obj_prev))) break;
        obj_prev = obj;
    }

    return Rcpp::List::create(
        Rcpp::Named("status")     = 3,
        Rcpp::Named("iterations") = iter,
        Rcpp::Named("M") = M,
        Rcpp::Named("S") = arma::exp(logS)
    );
}
