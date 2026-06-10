#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "nlopt_impl.h"

// ---------------------------------------------------------------------------------------
// Spherical covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_spherical(
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
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S);

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    const double w_bar = accu(w);
    const arma::uword p = Y.n_cols;
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat Xw  = X.each_col() % w;
    const arma::mat P_X = (X.n_cols > 0) ? arma::solve(X.t() * Xw, Xw.t()) : arma::mat(0, Y.n_rows);

    // E-step: M_full is the NLOPT parameter; B and sigma2 profiled at each eval
    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M_full = metadata.map<M_ID>(par);
        const arma::mat logS2  = metadata.map<S_ID>(par);
        const arma::mat S2     = arma::exp(logS2);
        const arma::mat B      = P_X * M_full;
        const arma::mat M_res  = M_full - X * B;
        const double sigma2    = accu(arma::diagmat(w) * (arma::pow(M_res, 2) + S2)) / (double(p) * w_bar);
        arma::mat gM, gS;
        const double obj = spherical_cov_obj_grad_impl(M_res, O + M_full, S2, logS2,
                                                        1./sigma2, 0.5 * double(p) * w_bar * std::log(sigma2),
                                                        Y, w, gM, gS);
        metadata.map<M_ID>(grad) = gM;
        metadata.map<S_ID>(grad) = gS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M      = metadata.copy<M_ID>(parameters.data());  // M_full
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat S      = arma::exp(0.5 * logS2);
    arma::mat B      = P_X * M;
    arma::mat M_res  = M - X * B;
    const double sigma2 = accu(diagmat(w) * (pow(M_res, 2) + S2)) / (double(p) * w_bar);
    arma::sp_mat Sigma(p, p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p, p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    arma::mat Z = O + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M_res, 2) + S2) / sigma2 + 0.5 * (logS2 - log(sigma2)), 1) + ki(Y);

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
            Rcpp::Named("backend", "nlopt"),
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

    const arma::mat XB    = X * B;  // B is fixed; precompute XB for M_res = M - XB
    const double omega2   = Omega(0, 0);  // fixed precision = 1/sigma2

    // Vestep: M_full is the NLOPT parameter; B and Omega fixed by the caller
    auto objective_and_grad = [&](const double * params, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);
        const arma::mat S2    = arma::exp(logS2);
        const arma::mat M_res = M - XB;
        const double penalty  = 0.5 * omega2 * accu(arma::diagmat(w) * (arma::pow(M_res, 2) + S2));
        arma::mat gM, gS;
        const double obj = spherical_cov_obj_grad_impl(M_res, O + M, S2, logS2,
                                                        omega2, penalty, Y, w, gM, gS);
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
    arma::mat loglik = sum(Y % Z - A - 0.5 * (pow(M_res, 2) + S2) * omega2 + 0.5 * (logS2 + log(omega2)), 1) + ki(Y);

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
