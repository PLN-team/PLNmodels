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
arma::mat optim_zipln_R_var(
    const arma::mat & Y,  // responses (n,p)
    const arma::mat & X,  // covariates (n,d)
    const arma::mat & O,  // offsets (n,p)
    const arma::mat & M,  // (n,p)
    const arma::mat & S2, // (n,p) variational variance
    const arma::mat & Pi, // (d,p)
    const arma::mat & B   // covariates (n,d)
) {
    arma::mat A = exp(O + M + 0.5 * S2);
    arma::mat R = 1. / (1. + exp(-(A + logit(Pi))));
    R %= arma::conv_to<arma::mat>::from(Y < 0.5);
    return R;
}

double phi (double mu, double sigma2) {
  double W = lambertW0_CS(sigma2 * exp(mu)) ;
  return(exp(-(pow(W, 2) + 2 * W) / (2 * sigma2)) / sqrt(1 + W)) ;
}

// [[Rcpp::export]]
arma::mat optim_zipln_R_exact (
    const arma::mat & Y,  // covariates (n,d)
    const arma::mat & X,  // covariates (n,d)
    const arma::mat & O,  // offsets (n,p)
    const arma::mat & M,  // (n,p)
    const arma::mat & S2, // (n,p) variational variance
    const arma::mat & Pi, // (n,p)
    const arma::mat & B   // covariates (n,d)
) {

  arma::mat XB = X * B;
  arma::mat M_mu = M - XB;
  const int n = (int)M.n_rows;
  const int p = (int)M.n_cols;
  arma::vec diag_Sigma = (sum(M_mu % M_mu, 0) + sum(S2, 0)).t() / double(n);
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

// ---------------------------------------------------------------------------------------
// Joint VE step for (M, ψ=log(S²), R) — nlopt backend.
//
// R is fixed at its exact conditional optimum R* = σ(A₀ + logit(Pi)) [Y>0 → 0]
// computed once from the initial (M, S²) before the nlopt solve.  Fixing R
// during the inner solve keeps the (objective, gradient) pair consistent, which
// is required by nlopt's line-search.  The final R is recomputed from the
// optimised (M, S²) before returning.
//
// B is kept fixed at the value from the preceding M-step (= P_X * M_prev).
// Dynamic B-profiling (like the Newton case) was tested but has no effect here:
// at the start of each VE step B is already at its optimum, so H*M ≈ H*M_prev
// throughout the small nlopt steps — the profiling adds O(dnp) cost per evaluation
// with no improvement in loglik.
//
// Interface: takes Pi (ZI structural probability), returns M, S², R.

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_nlopt(
    const arma::mat & init_M,  // (n,p)
    const arma::mat & init_S2, // (n,p) variational variance
    const arma::mat & Y,       // responses (n,p)
    const arma::mat & X,       // covariates (n,d)
    const arma::mat & O,       // offsets (n,p)
    const arma::mat & Pi,      // (n,p) ZI structural probability
    const arma::mat & B,       // (d,p)
    const arma::mat & Omega,   // (p,p)
    const Rcpp::List & configuration
) {
    const arma::mat psi_init   = arma::log(init_S2);
    const auto metadata        = tuple_metadata(init_M, psi_init);
    enum { M_ID, PSI_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>  (parameters.data()) = init_M;
    metadata.map<PSI_ID>(parameters.data()) = psi_init;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    const arma::mat X_B        = X * B;
    const arma::vec diag_Omega = diagvec(Omega);
    const arma::mat logit_Pi   = logit(Pi);
    const arma::mat Y_zero     = arma::conv_to<arma::mat>::from(Y < 0.5);

    // R is fixed at its exact conditional optimum given the initial (M, S2).
    // This ensures a consistent (objective, gradient) pair for nlopt.
    // After convergence the final R is recomputed from the optimised (M, S2).
    const arma::mat A0      = arma::exp(O + init_M + 0.5 * init_S2);
    const arma::mat R0      = (1.0 / (1.0 + arma::exp(-(A0 + logit_Pi)))) % Y_zero;
    const arma::mat one_m_R = 1.0 - R0;

    auto objective_and_grad = [&metadata, &Y, &O, &X_B, &Omega, &diag_Omega,
                                &one_m_R](
            const double * params, double * grad) -> double {
        const arma::mat M   = metadata.map<M_ID>  (params);
        const arma::mat psi = metadata.map<PSI_ID>(params);
        const arma::mat S2  = arma::exp(psi);
        const arma::mat A   = arma::exp(O + M + 0.5 * S2);
        const arma::mat M_mu       = M - X_B;
        const arma::mat M_mu_Omega = M_mu * Omega;

        double objective = - accu(one_m_R % (Y % M - A))
                         + 0.5 * accu(M_mu_Omega % M_mu)
                         + 0.5 * dot(diag_Omega, sum(S2, 0))
                         - 0.5 * accu(psi);

        metadata.map<M_ID>  (grad) = M_mu_Omega + one_m_R % (A - Y);
        metadata.map<PSI_ID>(grad) = 0.5 * (S2.each_row() % diag_Omega.t() + one_m_R % S2 % A - 1.);

        return objective;
    };

    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    const arma::mat M    = metadata.copy<M_ID>  (parameters.data());
    const arma::mat psi  = metadata.copy<PSI_ID>(parameters.data());
    const arma::mat S2   = arma::exp(psi);
    const arma::mat A    = arma::exp(O + M + 0.5 * S2);
    const arma::mat Rfin = (1.0 / (1.0 + arma::exp(-(A + logit_Pi)))) % Y_zero;
    // B_new = (X'X)^{-1} X' M — profiled analytically from the converged M.
    // Returned so the caller can skip the separate B M-step.
    const arma::mat B_new = (X.n_cols > 0)
        ? arma::solve(X.t() * X, X.t() * M, arma::solve_opts::likely_sympd)
        : arma::mat(0, (arma::uword)M.n_cols);
    return Rcpp::List::create(
        Rcpp::Named("status")     = static_cast<int>(result.status),
        Rcpp::Named("iterations") = result.nb_iterations,
        Rcpp::Named("M")  = M,
        Rcpp::Named("S2") = S2,
        Rcpp::Named("R")  = Rfin,
        Rcpp::Named("B")  = B_new
    );
}

// ---------------------------------------------------------------------------------------
// Joint VE Newton optimizer for (M, ψ=log(S²), R) — no NLopt dependency.
//
// B-profiling (envelope over regression coefficients, same as PLN builtin):
//   B_opt(M) = (X'X)^{-1} X' M  (closed-form WLS at each step)
//   M_res = (I − H) M            (in null(X') throughout)
// Each Newton step is projected onto null(X') via Q_step = step_M − X*(P_X*step_M),
// so M_res stays in null(X') and B is implicitly updated at every iteration.
// Z = O + M_full is tracked separately from M_res to correctly update A.
// By the envelope theorem the gradient w.r.t. M_full equals the gradient in M_res coords.
//
// R is the exact conditional optimum given (M, ψ): R_ij* = σ(A_ij + logit(π_ij)) for Y=0.
// Since ∂f/∂R = 0 at R*, updating R at each Newton iteration does not modify the gradient
// formula for (M, ψ) — only (1-R) changes, tightening the VE step.
//
// Per Newton iteration:
//   0. R  ← σ(A + logit(Pi)), zeroed where Y > 0          [exact VE for R, O(np)]
//   1. 2×2 Newton step for (M_res, ψ) with cross-term H_{Mψ}  [joint step]
//   1b. Project Newton step onto null(X'): Q_step = step_M − X*(P_X*step_M)  [B-profiling]
//   2. Joint Armijo: Z (M_full) moves by step_M; M_res moves by Q_step       [line search]
//
// M_res·Ω cached and updated incrementally (avoids O(np²) per trial).

// [[Rcpp::export]]
Rcpp::List ve_step_zipln_newton(
    const arma::mat & init_M,    // (n,p) M_full = X*B + M_res
    const arma::mat & init_S2,   // (n,p) variational variance
    const arma::mat & Y,         // (n,p)
    const arma::mat & X,         // (n,d)
    const arma::mat & O,         // (n,p)
    const arma::mat & Pi,        // (n,p) ZI structural probability (model param, fixed)
    const arma::mat & B,         // (d,p) — used to initialise M_res; then profiled out
    const arma::mat & Omega,     // (p,p)
    const int    maxiter,
    const double ftol_rel
) {
    const arma::mat XB         = X * B;
    const arma::vec diag_Omega = arma::diagvec(Omega);
    const arma::mat omega_d    = arma::ones(init_M.n_rows, 1) * diag_Omega.t(); // (n,p)
    const arma::mat logit_Pi   = logit(Pi);
    const arma::mat Y_zero     = arma::conv_to<arma::mat>::from(Y < 0.5);  // 1.0 where Y==0

    // B-profiling: P_X = (X'X)^{-1} X' (d×n). Projects Newton steps onto null(X')
    // so that M_res = (I−H)*M stays in null(X') and B is implicitly re-profiled.
    // Since B (from the M-step) = P_X*init_M, M_res = (I−H)*init_M ∈ null(X').
    const bool      do_profile = (X.n_cols > 0);
    const arma::mat P_X = do_profile
        ? arma::solve(X.t() * X, X.t(), arma::solve_opts::likely_sympd)
        : arma::mat(0, (arma::uword)init_M.n_rows);

    arma::mat M_res = init_M - XB;         // (I−H)*init_M: in null(X') by construction
    arma::mat psi   = arma::log(init_S2);
    arma::mat MO    = M_res * Omega;
    arma::mat S2    = arma::exp(psi);
    arma::mat Z     = O + init_M;          // Z = O + M_full; tracked separately from M_res
    arma::mat A     = arma::trunc_exp(Z + 0.5 * S2);

    // R and (1-R): updated at each Newton iteration, frozen during line search.
    arma::mat R     = arma::zeros(arma::size(A));
    arma::mat one_m_R(arma::size(A));

    // Objective: data term uses Z = O + M_full (captures B update via Z);
    // KL term uses M_res (projected, in null(X')).
    // Y*O is a constant offset — it cancels in all Armijo and convergence comparisons.
    auto obj_fun = [&](const arma::mat & Z_, const arma::mat & Mres,
                       const arma::mat & MO_, const arma::mat & psi_) -> double {
        const arma::mat S2_ = arma::exp(psi_);
        const arma::mat A_  = arma::trunc_exp(Z_ + 0.5 * S2_);
        return - arma::accu(Y % Z_) + arma::accu(one_m_R % A_)
               + 0.5 * arma::accu(Mres % MO_)
               + 0.5 * arma::dot(diag_Omega, arma::sum(S2_, 0))
               - 0.5 * arma::accu(psi_);
    };

    double obj_prev = 0.0;
    int iter = 0;

    for (iter = 0; iter < maxiter; iter++) {
        // ----------------------------------------------------------------
        // Step 0 — exact R update: R* = σ(A + logit(Pi)), R[Y>0] = 0
        // ∂f/∂R = 0 at R*, so the gradient w.r.t. (M,ψ) is unaffected.
        // ----------------------------------------------------------------
        R       = 1.0 / (1.0 + arma::exp(-(A + logit_Pi)));
        R      %= Y_zero;
        one_m_R = 1.0 - R;
        obj_prev = obj_fun(Z, M_res, MO, psi);

        // ----------------------------------------------------------------
        // Step 1 — joint 2×2 Newton step for (M_res, ψ)
        // ----------------------------------------------------------------
        const arma::mat one_m_R_A = one_m_R % A;

        const arma::mat h_mm = one_m_R_A + omega_d;
        const arma::mat h_mp = 0.5 * S2 % one_m_R_A;
        const arma::mat h_pp = h_mp % (1.0 + 0.5 * S2) + 0.5 * S2 % omega_d;

        const arma::mat grad_M   = MO + one_m_R_A - Y;
        const arma::mat grad_psi = h_mp + 0.5 * (S2 % omega_d - 1.0);

        arma::mat det = h_mm % h_pp - h_mp % h_mp;
        det.clamp(1e-20, arma::datum::inf);

        const arma::mat step_M   = (h_pp % grad_M   - h_mp % grad_psi) / det;
        const arma::mat step_psi = (h_mm % grad_psi - h_mp % grad_M  ) / det;

        // ----------------------------------------------------------------
        // Step 1b — B-profiling: project Newton step onto null(X').
        // Q_step moves M_res (KL term); step_M moves Z = O + M_full (data term A).
        // By the envelope theorem the gradient w.r.t. M_full = gradient in M_res coords.
        // ----------------------------------------------------------------
        const arma::mat Q_step = do_profile ? step_M - X * (P_X * step_M) : step_M;
        const arma::mat dMO    = Q_step * Omega;  // incremental MO update

        // ----------------------------------------------------------------
        // Step 2 — joint Armijo: Z moves by step_M (for A); M_res by Q_step (for KL).
        // ----------------------------------------------------------------
        double slope = -arma::accu(grad_M % Q_step) - arma::accu(grad_psi % step_psi);
        if (slope >= 0.0)
            slope = -(arma::accu(arma::square(grad_M)) + arma::accu(arma::square(grad_psi)));

        constexpr double c1 = 1e-4;
        double alpha        = 1.0;
        for (int ls = 0; ls < 20; ++ls) {
            if (obj_fun(Z     - alpha * step_M,   // full step: tracks M_full for A
                        M_res - alpha * Q_step,    // projected: KL residual stays in null(X')
                        MO    - alpha * dMO,
                        psi   - alpha * step_psi) <= obj_prev + c1 * alpha * slope) break;
            alpha *= 0.5;
        }

        M_res -= alpha * Q_step;
        Z     -= alpha * step_M;     // Z = O + M_full_new; B implicitly = P_X*(Z-O)
        MO    -= alpha * dMO;
        psi   -= alpha * step_psi;
        S2     = arma::exp(psi);
        A      = arma::trunc_exp(Z + 0.5 * S2);

        // ----------------------------------------------------------------
        // Convergence: compare objective before and after this Newton step
        // (R is frozen during the step; obj_prev was set at top of this iter)
        // ----------------------------------------------------------------
        const double obj = obj_fun(Z, M_res, MO, psi);
        if (std::abs(obj - obj_prev) < ftol_rel * (1.0 + std::abs(obj_prev))) break;
    }

    const arma::mat M_full = Z - O;
    // B profiled from final M_full: reuse P_X already computed above.
    const arma::mat B_new = do_profile ? P_X * M_full : B;
    return Rcpp::List::create(
        Rcpp::Named("status")     = 3,
        Rcpp::Named("iterations") = iter,
        Rcpp::Named("M")          = M_full,
        Rcpp::Named("S2")         = S2,
        Rcpp::Named("R")          = R,
        Rcpp::Named("B")          = B_new
    );
}
