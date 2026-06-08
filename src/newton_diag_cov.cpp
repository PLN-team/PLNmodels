#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "newton_impl.h"

// ---------------------------------------------------------------------------------------
// Diagonal covariance PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List newton_optimize_diagonal(
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
    DiagonalCovTraits::State state(M, S2, w, w_bar);

    return newton_optimize_impl<DiagonalCovTraits>(Y, X, O, w, B, M, S, state, maxiter, ftol, max_em, em_tol);
}

// ---------------------------------------------------------------------------------------
// VE diagonal — coordinate-Newton (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List newton_optimize_vestep_diagonal(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p) fixed
    const arma::mat & Omega,   // (p,p) diagonal, fixed
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

    const int n = Y.n_rows;
    const double c1 = 1e-4;
    const arma::rowvec omega2 = arma::diagvec(Omega).t();
    const arma::mat ones_row  = arma::ones(n, 1);
    const arma::mat XB        = X * B;

    arma::mat logS = arma::log(S);
    arma::mat S2   = S % S;

    std::vector<double> objective_vec;
    double obj_prev   = arma::datum::inf;
    int    total_iter = 0;

    for (int it = 0; it < maxiter; it++) {
        S2 = S % S;
        arma::mat Z = O + XB + M;
        arma::mat A = arma::exp(Z + 0.5 * S2);

        // ---- Diagonal Newton step for M ----
        arma::mat grad_M = M.each_row() % omega2 + A - Y; grad_M.each_col() %= w;
        arma::mat hess_M = ones_row * omega2 + A;          hess_M.each_col() %= w;
        hess_M.clamp(1e-10, arma::datum::inf);
        arma::mat step_M = grad_M / hess_M;
        double f0_M    = arma::accu(w.t() * (A - Y % Z))
                       + 0.5 * arma::as_scalar((w.t() * (M % M)) * omega2.t());
        double slope_M = -arma::accu(grad_M % step_M);
        double alpha_M = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            arma::mat Mt = M - alpha_M * step_M;
            arma::mat Zt = Z - alpha_M * step_M;
            arma::mat At = arma::exp(Zt + 0.5 * S2);
            if (arma::accu(w.t() * (At - Y % Zt))
                + 0.5 * arma::as_scalar((w.t() * (Mt % Mt)) * omega2.t())
                <= f0_M + c1 * alpha_M * slope_M) break;
            alpha_M *= 0.5;
        }
        M -= alpha_M * step_M;
        Z  = O + XB + M;

        // ---- Fixed-point update for logS ----
        fixed_point_logS(logS, S, S2, Z, A, ones_row * omega2);

        // ---- Objective for convergence ----
        A = arma::exp(Z + 0.5 * S2);
        double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * arma::trunc_log(S2)))
                   + 0.5 * arma::as_scalar((w.t() * (M % M + S2)) * omega2.t());
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
    }

    // ---- Final output ----
    S2 = S % S;
    arma::mat Z = O + XB + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec omega2_v = omega2.t();
    arma::vec loglik = arma::sum(Y % Z - A + 0.5 * arma::log(S2), 1)
                     - 0.5 * (M % M + S2) * omega2_v
                     + 0.5 * arma::accu(arma::log(omega2_v)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S")  = S,
        Rcpp::Named("Ji") = Ji,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     3            ),
            Rcpp::Named("backend",    "newton"     ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}
