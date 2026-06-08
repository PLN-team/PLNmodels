#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "newton_impl.h"

// ---------------------------------------------------------------------------------------
// Diagonal covariance

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_diagonal(
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

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));
    const double w_bar = accu(w);

    const arma::mat Xw = X.each_col() % w;   // fixed: precomputed once

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Xw, &Y, &w, &w_bar, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::rowvec diag_sigma = w.t() * (M % M + S2) / w_bar;
        double objective = accu(diagmat(w) * (A - Y % Z - 0.5 * log(S2))) + 0.5 * w_bar * accu(log(diag_sigma));

        metadata.map<B_ID>(grad) = Xw.t() * (A - Y);
        metadata.map<M_ID>(grad) = diagmat(w) * ((M.each_row() / diag_sigma) + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % pow(diag_sigma, -1) + S % A - pow(S, -1));

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
    arma::rowvec sigma2 = w.t() * (M % M + S2) / w_bar;
    arma::vec omega2 = pow(sigma2.t(), -1);
    arma::sp_mat Sigma(Y.n_cols, Y.n_cols);
    Sigma.diag() = sigma2;
    arma::sp_mat Omega(Y.n_cols, Y.n_cols);
    Omega.diag() = omega2;
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik =
        sum(Y % Z - A + 0.5 * log(S2), 1) - 0.5 * (pow(M, 2) + S2) * omega2 + 0.5 * sum(log(omega2)) + ki(Y);

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
// Diagonal covariance PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_newton_diagonal(
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
Rcpp::List nlopt_optimize_vestep_newton_diagonal(
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

// ---------------------------------------------------------------------------------------
// VE diagonal

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_diagonal(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p)
    const arma::mat & Omega,   // (p,p)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
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
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB_diag  = X * B;
    const arma::vec omega2_v = arma::diagvec(Omega);  // fixed: Omega not optimized in vestep

    // Optimize
    auto objective_and_grad = [&metadata, &O, &XB_diag, &Y, &w, &omega2_v, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + XB_diag + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double objective =
            accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * as_scalar(w.t() * (pow(M, 2) + S2) * omega2_v) ;

        metadata.map<M_ID>(grad) = diagmat(w) * (M * arma::diagmat(omega2_v) + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % omega2_v.t() + S % A - pow(S, -1));

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::vec omega2 = Omega.diag();
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik =
      sum(Y % Z - A + 0.5 * log(S2), 1) - 0.5 * (pow(M, 2) + S2) * omega2 + 0.5 * sum(log(omega2)) + ki(Y);

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
