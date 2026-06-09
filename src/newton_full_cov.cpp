#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"
#include "newton_impl_alt.h"

// ---------------------------------------------------------------------------------------
// Full covariance PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List newton_optimize_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter = config.containsElementNamed("maxeval")     ? Rcpp::as<int>(config["maxeval"])        : 200;
    const double ftol    = config.containsElementNamed("ftol_rel")    ? Rcpp::as<double>(config["ftol_rel"])    : 1e-8;
    const int    max_em  = config.containsElementNamed("max_em_iter") ? Rcpp::as<int>(config["max_em_iter"])    : 50;
    const double em_tol  = config.containsElementNamed("em_ftol")     ? Rcpp::as<double>(config["em_ftol"])     : 1e-8;

    const double w_bar = arma::accu(w);
    arma::mat S2 = S % S;
    FullCovTraits::State state(M, S2, w, w_bar);

    return newton_optimize_impl<FullCovTraits>(Y, X, O, w, B, M, S, state, maxiter, ftol, max_em, em_tol);
}

// ---------------------------------------------------------------------------------------
// Full covariance PLN — alternative EM: closed-form B in M-step, no newton_step_B in inner loop

// [[Rcpp::export]]
Rcpp::List newton_optimize_full_alt(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter = config.containsElementNamed("maxeval")     ? Rcpp::as<int>(config["maxeval"])        : 200;
    const double ftol    = config.containsElementNamed("ftol_rel")    ? Rcpp::as<double>(config["ftol_rel"])    : 1e-8;
    const int    max_em  = config.containsElementNamed("max_em_iter") ? Rcpp::as<int>(config["max_em_iter"])    : 50;
    const double em_tol  = config.containsElementNamed("em_ftol")     ? Rcpp::as<double>(config["em_ftol"])     : 1e-8;

    const double w_bar = arma::accu(w);
    arma::mat S2 = S % S;
    FullCovTraits::State state(M, S2, w, w_bar);

    return newton_optimize_alt_impl<FullCovTraits>(Y, X, O, w, B, M, S, state, maxiter, ftol, max_em, em_tol);
}

// ---------------------------------------------------------------------------------------
// VE full covariance — coordinate-Newton (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List newton_optimize_vestep_full(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p) fixed
    const arma::mat & Omega,   // (p,p) fixed
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 200;
    const double ftol    = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-8;

    FullCovTraits::State state(Omega);
    return newton_vestep_impl<FullCovTraits>(Y, X, O, w, M, S, B, state, maxiter, ftol);
}
