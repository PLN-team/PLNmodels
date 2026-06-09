#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Spherical covariance PLN — nlopt/CCSAQ optimizer

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_spherical(
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

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    const double w_bar = accu(w);
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat Xw = X.each_col() % w;   // fixed: precomputed once

    auto objective_and_grad = [&metadata, &O, &X, &Xw, &Y, &w, &w_bar, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat B     = metadata.map<B_ID>(params);
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);

        arma::mat S2 = arma::exp(logS2);
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
        // -½ log(S²) → -½ logS2
        double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2)) + 0.5 * (double(p) * w_bar) * log(sigma2) ;

        metadata.map<B_ID>(grad) = Xw.t() * (A - Y);
        metadata.map<M_ID>(grad) = diagmat(w) * (M / sigma2 + A - Y);
        // grad_logS2 = ½ w ⊙ (S²/σ² + S²⊙A − 1)
        metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2 / sigma2 + S2 % A - 1.);

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Variational parameters
    arma::mat M     = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2    = arma::exp(logS2);
    arma::mat S     = arma::exp(0.5 * logS2);
    // Regression parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    // Variance parameters
    const arma::uword p = Y.n_cols;
    const double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
    arma::sp_mat Sigma(p,p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p,p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    // Element-wise log-likelihood   [log(S²/σ²) = logS2 - log(σ²)]
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2) / sigma2 + 0.5 * (logS2 - log(sigma2)), 1) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Ji", Ji),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("objective", objective_vec),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
      );
}

// ---------------------------------------------------------------------------------------
// Spherical covariance PLN — profiled-B nlopt: B removed from parameter vector, closed-form per eval

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_spherical_alt(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    const auto init_B = Rcpp::as<arma::mat>(params["B"]);
    const auto init_M = Rcpp::as<arma::mat>(params["M"]);
    const auto init_S = Rcpp::as<arma::mat>(params["S"]);

    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = X * init_B + init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S);

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    const double w_bar = accu(w);
    const arma::uword p = Y.n_cols;
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat Xw  = X.each_col() % w;
    const arma::mat P_X = arma::solve(X.t() * Xw, Xw.t());

    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M_full = metadata.map<M_ID>(par);
        const arma::mat logS2  = metadata.map<S_ID>(par);
        arma::mat S2    = arma::exp(logS2);
        arma::mat B     = P_X * M_full;
        arma::mat M_res = M_full - X * B;
        arma::mat Z     = O + M_full;
        arma::mat A     = exp(Z + 0.5 * S2);
        double sigma2   = accu(diagmat(w) * (pow(M_res, 2) + S2)) / (double(p) * w_bar);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2)) + 0.5 * double(p) * w_bar * log(sigma2);
        // gradient for M_full = gradient for M_res (envelope theorem for B and sigma2)
        metadata.map<M_ID>(grad) = diagmat(w) * (M_res / sigma2 + A - Y);
        metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2 / sigma2 + S2 % A - 1.);
        objective_vec.push_back(objective);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M_full = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat S      = arma::exp(0.5 * logS2);
    arma::mat B      = P_X * M_full;
    arma::mat M      = M_full - X * B;
    const double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar);
    arma::sp_mat Sigma(p, p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p, p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    arma::mat Z = O + M_full;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2) / sigma2 + 0.5 * (logS2 - log(sigma2)), 1) + ki(Y);

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
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt_alt"),
            Rcpp::Named("objective", objective_vec),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE spherical — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_spherical(
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
    const auto init_S = Rcpp::as<arma::mat>(params["S"]); // (n)

    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S); // pack logS2

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB_sph = X * B;   // fixed: B not optimized in vestep

    auto objective_and_grad = [&metadata, &O, &XB_sph, &Y, &w, &Omega, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);

        arma::mat S2 = arma::exp(logS2);
        arma::mat Z = O + XB_sph + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double n_sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) ;
        double omega2 = Omega(0, 0);
        // -½ log(S²) → -½ logS2
        double objective = accu(w.t() * (A - Y % Z - 0.5 * logS2)) + 0.5 * n_sigma2 * omega2;

        metadata.map<M_ID>(grad) = diagmat(w) * (M / omega2 + A - Y);
        // grad_logS2 = ½ w ⊙ (S²/ω² + S²⊙A − 1)
        metadata.map<S_ID>(grad) = 0.5 * diagmat(w) * (S2 / omega2 + S2 % A - 1.);

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M     = metadata.copy<M_ID>(parameters.data());
    arma::mat logS2 = metadata.copy<S_ID>(parameters.data());
    arma::mat S2    = arma::exp(logS2);
    arma::mat S     = arma::exp(0.5 * logS2);
    double omega2 = Omega(0, 0);
    // Element-wise log-likelihood  [log(S²·ω²) = logS2 + log(ω²)]
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2) * omega2 + 0.5 * (logS2 + log(omega2)), 1) + ki(Y);

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
