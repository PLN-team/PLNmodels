#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Shared inner computation for diagonal-covariance NLOPT objective and gradients.
// Z = O + X*B + M_res and S2 = exp(logS2) must be pre-computed by the caller.
// inv_sigma2: row vector of element-wise precisions (profiled 1/diag_sigma or fixed omega2).
// penalty: KL variance term (log-det for profiled E-step; quadratic for fixed-omega vestep).
static double diag_cov_obj_grad_impl(
    const arma::mat & M_res, const arma::mat & Z,
    const arma::mat & S2,   const arma::mat & logS2,
    const arma::rowvec & inv_sigma2, double penalty,
    const arma::mat & Y, const arma::vec & w,
    arma::mat & grad_M, arma::mat & grad_S
) {
    const arma::mat A = arma::exp(Z + 0.5 * S2);
    grad_M = arma::diagmat(w) * (M_res.each_row() % inv_sigma2 + A - Y);
    grad_S = 0.5 * arma::diagmat(w) * (S2.each_row() % inv_sigma2 + S2 % A - 1.);
    return accu(w.t() * (A - Y % Z - 0.5 * logS2)) + penalty;
}

// ---------------------------------------------------------------------------------------
// Diagonal covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_diagonal(
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
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));
    const double w_bar = accu(w);

    const arma::mat Xw  = X.each_col() % w;
    const arma::mat P_X = arma::solve(X.t() * Xw, Xw.t());

    // E-step: M_full is the NLOPT parameter; B and diag_sigma profiled at each eval
    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M_full    = metadata.map<M_ID>(par);
        const arma::mat logS2     = metadata.map<S_ID>(par);
        const arma::mat S2        = arma::exp(logS2);
        const arma::mat B         = P_X * M_full;
        const arma::mat M_res     = M_full - X * B;
        const arma::rowvec diag_sigma = w.t() * (M_res % M_res + S2) / w_bar;
        const arma::rowvec inv_sigma2 = arma::pow(diag_sigma, -1);
        arma::mat gM, gS;
        const double obj = diag_cov_obj_grad_impl(M_res, O + M_full, S2, logS2,
                                                   inv_sigma2, 0.5 * w_bar * accu(arma::log(diag_sigma)),
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
    arma::rowvec sigma2 = w.t() * (M_res % M_res + S2) / w_bar;
    arma::vec omega2    = pow(sigma2.t(), -1);
    arma::sp_mat Sigma(Y.n_cols, Y.n_cols); Sigma.diag() = sigma2.t();
    arma::sp_mat Omega(Y.n_cols, Y.n_cols); Omega.diag() = omega2;
    arma::mat Z = O + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A + 0.5 * logS2, 1) - 0.5 * (pow(M_res, 2) + S2) * omega2
                     + 0.5 * sum(log(omega2)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("M", M),        // M_full
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
// VE diagonal — nlopt/CCSAQ (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_diagonal(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p)
    const arma::mat & Omega,   // (p,p)
    const Rcpp::List & config  // List of config values
) {
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

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB           = X * B;  // B is fixed; precompute XB for M_res = M - XB
    const arma::rowvec omega2_v  = arma::diagvec(Omega).t();  // fixed precision, as row vector

    // Vestep: M_full is the NLOPT parameter; B and Omega fixed by the caller
    auto objective_and_grad = [&](const double * params, double * grad) -> double {
        const arma::mat M     = metadata.map<M_ID>(params);
        const arma::mat logS2 = metadata.map<S_ID>(params);
        const arma::mat S2    = arma::exp(logS2);
        const arma::mat M_res = M - XB;
        const double penalty  = 0.5 * as_scalar(w.t() * (arma::pow(M_res, 2) + S2) * omega2_v.t());
        arma::mat gM, gS;
        const double obj = diag_cov_obj_grad_impl(M_res, O + M, S2, logS2, omega2_v, penalty, Y, w, gM, gS);
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
    arma::vec omega2 = Omega.diag();
    arma::mat loglik =
      sum(Y % Z - A + 0.5 * logS2, 1) - 0.5 * (pow(M_res, 2) + S2) * omega2 + 0.5 * sum(log(omega2)) + ki(Y);

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
