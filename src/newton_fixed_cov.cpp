#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"

// ---------------------------------------------------------------------------------------
// Fixed inverse covariance (Omega provided externally) PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List newton_optimize_fixed(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S, Omega)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B     = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M     = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S     = Rcpp::as<arma::mat>(params["S"]);
    arma::mat Omega = Rcpp::as<arma::mat>(params["Omega"]);

    const int    maxiter = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 200;
    const double ftol    = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-8;

    arma::mat S2 = S % S;
    FixedCovTraits::State state(Omega);

    // Fixed covariance has no EM loop (has_em = false); pass max_em=1, em_tol=0
    return newton_optimize_impl<FixedCovTraits>(Y, X, O, w, B, M, S, state, maxiter, ftol, 1, 0.0);
}
