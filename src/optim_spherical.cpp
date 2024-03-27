#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Spherical covariance

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
    metadata.map<S_ID>(parameters.data()) = init_S;

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    const double w_bar = accu(w);
    std::vector<double> objective_vec ;

    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &w_bar, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * (double(p) * w_bar) * log(sigma2) ;

        metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);
        metadata.map<M_ID>(grad) = diagmat(w) * (M / sigma2 + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S / sigma2 + S % A - pow(S, -1)) ;

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    // Variance parameters
    const arma::uword p = Y.n_cols;
    const double sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) / (double(p) * w_bar) ;
    arma::sp_mat Sigma(p,p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p,p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2 ) / sigma2 + 0.5 * log(S2 / sigma2), 1) + ki(Y);

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
// VE spherical

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
    metadata.map<S_ID>(parameters.data()) = init_S;


    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;

    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &B, &Omega, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double n_sigma2 = accu(diagmat(w) * (pow(M, 2) + S2)) ;
        double omega2 = Omega(0, 0);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * n_sigma2 * omega2;

        metadata.map<M_ID>(grad) = diagmat(w) * (M / omega2 + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S / omega2 + S % A - pow(S, -1));

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    double omega2 = Omega(0, 0);
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M, 2) + S2 ) * omega2 + 0.5 * log(S2 * omega2), 1) + ki(Y);

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
