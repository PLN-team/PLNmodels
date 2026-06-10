#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "nlopt_impl.h"

// ---------------------------------------------------------------------------------------
// Full covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    const auto init_B = Rcpp::as<arma::mat>(params["B"]);
    const auto init_M = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S = Rcpp::as<arma::mat>(params["S"]);

    // Parameters: (M_full, logS2) — B is profiled out via closed form
    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S);

    const double w_bar = accu(w);
    const NewtonConfig cfg(config);

    // P_X = (X'WX)^{-1} X'W : d×n, precomputed once; B = P_X * M_full at each eval
    const arma::mat Xw  = X.each_col() % w;
    const arma::mat P_X = (X.n_cols > 0) ? arma::solve(X.t() * Xw, Xw.t()) : arma::mat(0, Y.n_rows);

    // Initial Omega: M_res = M_full - X*B
    arma::mat Omega;
    {
        const arma::mat S2_init   = init_S % init_S;
        const arma::mat M_res_init = init_M - X * init_B;
        arma::mat Sigma_init = (1./w_bar) * (M_res_init.t() * (M_res_init.each_col() % w) + diagmat(w.t() * S2_init));
        Omega = inv_sympd(Sigma_init);
    }

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iterations = 0;
    int last_status = 0;

    for (int em_iter = 0; em_iter < std::max(1, cfg.max_em); em_iter++) {
        auto optimizer = new_nlopt_optimizer(config, parameters.size());
        objective_vec.reserve(objective_vec.size() + nlopt_get_maxeval(optimizer.get()));
        const arma::vec Omega_diag = diagvec(Omega);

        // E-step: M_full is the NLOPT parameter; B profiled at each eval (envelope theorem)
        auto objective_and_grad = [&](const double * par, double * grad) -> double {
            const arma::mat M_full = metadata.map<M_ID>(par);
            const arma::mat logS2  = metadata.map<S_ID>(par);
            const arma::mat B      = P_X * M_full;
            const arma::mat M_res  = M_full - X * B;
            arma::mat gM, gS;
            const double obj = full_cov_obj_grad_impl(M_res, O + M_full, logS2, Omega, Omega_diag, Y, w, gM, gS);
            metadata.map<M_ID>(grad) = gM;
            metadata.map<S_ID>(grad) = gS;
            objective_vec.push_back(obj);
            return obj;
        };

        OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
        total_iterations += result.nb_iterations;
        last_status = static_cast<int>(result.status);

        arma::mat M_full = metadata.copy<M_ID>(parameters.data());
        arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
        arma::mat S2     = arma::exp(logS2);
        arma::mat B      = P_X * M_full;
        arma::mat M_res  = M_full - X * B;
        arma::mat Sigma  = (1./w_bar) * (M_res.t() * (M_res.each_col() % w) + diagmat(w.t() * S2));
        Omega = inv_sympd(Sigma);

        arma::mat Z = O + M_full;
        arma::mat A = exp(Z + 0.5 * S2);
        double elbo = accu(w.t() * (Y % Z - A + 0.5 * logS2))
                    - 0.5 * w_bar * real(log_det(Sigma));
        if (em_iter > 0 && converged(elbo, elbo_prev, cfg.em_tol)) break;
        elbo_prev = elbo;
    }

    arma::mat M      = metadata.copy<M_ID>(parameters.data());  // M_full
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat S      = arma::exp(0.5 * logS2);
    arma::mat B      = P_X * M;
    arma::mat M_res  = M - X * B;
    arma::mat Sigma  = (1./w_bar) * (M_res.t() * (M_res.each_col() % w) + diagmat(w.t() * S2));
    arma::mat Z      = O + M;
    arma::mat A      = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * logS2 - 0.5 * ((M_res * Omega) % M_res + S2 * diagmat(Omega)), 1)
                     + 0.5 * real(log_det(Omega)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji", Ji),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", last_status),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("objective", objective_vec),
            Rcpp::Named("iterations", total_iterations)
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE full covariance — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p)
    const arma::mat & Omega,   // (p,p)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(params["S"]); // (n,p)

    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S); // pack logS2

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB        = X * B;  // B is fixed; precompute XB for M_res = M - XB
    const arma::vec Omega_diag = diagvec(Omega);

    // Vestep: M_full is the NLOPT parameter; B and Omega fixed by the caller
    auto objective_and_grad = [&](const double * params, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);
        arma::mat gM, gS;
        const double obj = full_cov_obj_grad_impl(M - XB, O + M, logS2, Omega, Omega_diag, Y, w, gM, gS);
        metadata.map<M_ID>(grad) = gM;
        metadata.map<S_ID>(grad) = gS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M     = metadata.copy<M_ID>(parameters.data());  // M_full
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2    = arma::exp(logS2);
    arma::mat S     = arma::exp(0.5 * logS2);
    arma::mat M_res = M - XB;
    // Element-wise log-likelihood
    arma::mat Z = O + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * logS2 - 0.5 * ((M_res * Omega) % M_res + S2 * diagmat(Omega)), 1) +
      0.5 * real(log_det(Omega)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
      Rcpp::Named("M") = M,
      Rcpp::Named("S") = S,
      Rcpp::Named("Ji") = Ji,
      Rcpp::Named("monitoring", Rcpp::List::create(
          Rcpp::Named("status", static_cast<int>(result.status)),
          Rcpp::Named("backend", "nlopt"),
          Rcpp::Named("objective", objective_vec),
          Rcpp::Named("iterations", result.nb_iterations)
      ))
    );
}
