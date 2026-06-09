#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Full covariance PLN — nlopt/CCSAQ optimizer

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (d,p)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(params["S"]); // (n,p)

    const auto metadata = tuple_metadata(init_B, init_M, init_S);
    enum { B_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S); // pack logS2

    const double w_bar = accu(w);
    // VEM config — sensible defaults if not supplied by R
    const int    max_em_iter = config.containsElementNamed("max_em_iter") ? Rcpp::as<int>(config["max_em_iter"])       : 50;
    const double em_ftol     = config.containsElementNamed("em_ftol")     ? Rcpp::as<double>(config["em_ftol"])        : 1e-8;

    // Initial Omega: closed-form from initial M, S (one inv_sympd, outside the loop)
    arma::mat Omega;
    {
        const arma::mat S2_init = init_S % init_S;
        arma::mat Sigma_init = (1. / w_bar) * (init_M.t() * (init_M.each_col() % w) + diagmat(w.t() * S2_init));
        Omega = inv_sympd(Sigma_init);
    }

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iterations = 0;
    int last_status = 0;
    const arma::mat Xw = X.each_col() % w;   // fixed: precomputed once for all EM iterations

    for (int em_iter = 0; em_iter < std::max(1, max_em_iter); em_iter++) {
        // E-step: optimize B, M, S with Omega fixed — no inv_sympd inside gradient
        auto optimizer = new_nlopt_optimizer(config, parameters.size());
        objective_vec.reserve(objective_vec.size() + nlopt_get_maxeval(optimizer.get()));
        const arma::vec Omega_diag = diagvec(Omega);  // fixed per EM iteration

        auto objective_and_grad = [&metadata, &Y, &X, &Xw, &O, &w, &Omega, &Omega_diag, &objective_vec](const double * params, double * grad) -> double {
            const arma::mat B    = metadata.map<B_ID>(params);
            const arma::mat M    = metadata.map<M_ID>(params);
            const arma::mat logS2 = metadata.map<S_ID>(params); // logS2 = log(S²)
            arma::mat S2 = arma::exp(logS2);
            arma::mat Z = O + X * B + M;
            arma::mat A = exp(Z + 0.5 * S2);
            arma::mat MO = M * Omega;           // cached: reused in objective and M-gradient
            const arma::rowvec wS2 = w.t() * S2;

            // -½ log(S²) → -½ logS2 (no log call needed)
            double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2))
                             + 0.5 * (accu(MO % (M.each_col() % w)) + dot(Omega_diag, wS2.t()));

            metadata.map<B_ID>(grad) = Xw.t() * (A - Y);
            metadata.map<M_ID>(grad) = diagmat(w) * (MO + A - Y);
            // grad_logS2 = ½ w ⊙ (S²⊙(Ω_diag + A) − 1)
            metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2.each_row() % Omega_diag.t() + S2 % A - 1.);

            objective_vec.push_back(objective);
            return objective;
        };

        OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
        total_iterations += result.nb_iterations;
        last_status = static_cast<int>(result.status);

        // M-step: update Omega analytically (one inv_sympd per EM iteration)
        arma::mat M    = metadata.copy<M_ID>(parameters.data());
        arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
        arma::mat S2   = arma::exp(logS2);
        arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
        Omega = inv_sympd(Sigma);

        // ELBO after M-step: trace(Omega*nSigma) = w_bar*p, log_det(Omega) = -log_det(Sigma)
        arma::mat B = metadata.copy<B_ID>(parameters.data());
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double elbo = accu(w.t() * (Y % Z - A + 0.5 * logS2))
                    - 0.5 * w_bar * real(log_det(Sigma));

        if (em_iter > 0 && converged(elbo, elbo_prev, em_ftol)) break;
        elbo_prev = elbo;
    }

    // Final extraction
    arma::mat M    = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2   = arma::exp(logS2);
    arma::mat S    = arma::exp(0.5 * logS2);
    arma::mat B    = metadata.copy<B_ID>(parameters.data());
    arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    // Omega already updated from the last M-step
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * logS2 - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) +
                       0.5 * real(log_det(Omega)) + ki(Y);

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
// Full covariance PLN — profiled-B nlopt: B removed from parameter vector, closed-form per eval

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_full_alt(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    const auto init_B = Rcpp::as<arma::mat>(params["B"]);
    const auto init_M = Rcpp::as<arma::mat>(params["M"]);  // residual format
    const auto init_S = Rcpp::as<arma::mat>(params["S"]);

    // Parameters: (M_full, logS2) — B is profiled out via closed form
    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = X * init_B + init_M;  // M_full = XB + M_res
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S);

    const double w_bar    = accu(w);
    const int max_em_iter = config.containsElementNamed("max_em_iter") ? Rcpp::as<int>(config["max_em_iter"])    : 50;
    const double em_ftol  = config.containsElementNamed("em_ftol")     ? Rcpp::as<double>(config["em_ftol"])     : 1e-8;

    // P_X = (X'WX)^{-1} X'W : d×n, precomputed once; B = P_X * M_full at each eval
    const arma::mat Xw  = X.each_col() % w;
    const arma::mat P_X = arma::solve(X.t() * Xw, Xw.t());

    // Initial Omega from residual init_M
    arma::mat Omega;
    {
        const arma::mat S2_init = init_S % init_S;
        arma::mat Sigma_init = (1./w_bar) * (init_M.t() * (init_M.each_col() % w) + diagmat(w.t() * S2_init));
        Omega = inv_sympd(Sigma_init);
    }

    std::vector<double> objective_vec;
    double elbo_prev = -arma::datum::inf;
    int total_iterations = 0;
    int last_status = 0;

    for (int em_iter = 0; em_iter < std::max(1, max_em_iter); em_iter++) {
        auto optimizer = new_nlopt_optimizer(config, parameters.size());
        objective_vec.reserve(objective_vec.size() + nlopt_get_maxeval(optimizer.get()));
        const arma::vec Omega_diag = diagvec(Omega);

        auto objective_and_grad = [&](const double * par, double * grad) -> double {
            const arma::mat M_full = metadata.map<M_ID>(par);
            const arma::mat logS2  = metadata.map<S_ID>(par);
            arma::mat S2    = arma::exp(logS2);
            arma::mat B     = P_X * M_full;
            arma::mat M_res = M_full - X * B;
            arma::mat Z     = O + M_full;
            arma::mat A     = exp(Z + 0.5 * S2);
            arma::mat MO    = M_res * Omega;
            const arma::rowvec wS2 = w.t() * S2;
            double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2))
                             + 0.5 * (accu(MO % (M_res.each_col() % w)) + dot(Omega_diag, wS2.t()));
            // gradient for M_full = gradient for M_res (envelope theorem for B)
            metadata.map<M_ID>(grad) = diagmat(w) * (MO + A - Y);
            metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2.each_row() % Omega_diag.t() + S2 % A - 1.);
            objective_vec.push_back(objective);
            return objective;
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
        if (em_iter > 0 && converged(elbo, elbo_prev, em_ftol)) break;
        elbo_prev = elbo;
    }

    arma::mat M_full = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat S      = arma::exp(0.5 * logS2);
    arma::mat B      = P_X * M_full;
    arma::mat M      = M_full - X * B;  // M_res — residual format
    arma::mat Sigma  = (1./w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    arma::mat Z      = O + M_full;
    arma::mat A      = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * logS2 - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1)
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
            Rcpp::Named("backend", "nlopt_alt"),
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

    const arma::mat XB_vestep   = X * B;
    const arma::vec Omega_diag_v = diagvec(Omega);

    auto objective_and_grad = [&metadata, &O, &XB_vestep, &Y, &w, &Omega, &Omega_diag_v, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);

        arma::mat S2 = arma::exp(logS2);
        arma::mat Z = O + XB_vestep + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat MO = M * Omega;
        const arma::rowvec wS2 = w.t() * S2;
        double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2))
                         + 0.5 * (accu(MO % (M.each_col() % w)) + dot(Omega_diag_v, wS2.t()));

        metadata.map<M_ID>(grad) = diagmat(w) * (MO + A - Y);
        metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2.each_row() % Omega_diag_v.t() + S2 % A - 1.);

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M     = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2    = arma::exp(logS2);
    arma::mat S     = arma::exp(0.5 * logS2);
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * logS2 - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) +
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
