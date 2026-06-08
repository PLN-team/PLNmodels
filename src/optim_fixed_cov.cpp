#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Fixed inverse covariance (Omega)

// ---------------------------------------------------------------------------------------
// Fixed covariance PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_newton_fixed(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S, Omega)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y     = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X     = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O     = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w     = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B             = Rcpp::as<arma::mat>(params["B"]);
    arma::mat M             = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S             = Rcpp::as<arma::mat>(params["S"]);
    const arma::mat   Omega = Rcpp::as<arma::mat>(params["Omega"]);

    const int maxiter  = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 200;
    const double ftol  = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-8;

    const int n = Y.n_rows;
    const double w_bar = arma::accu(w);
    const double c1 = 1e-4;

    const arma::mat Xw  = X.each_col() % w;
    arma::mat Xw2 = X % X; Xw2.each_col() %= w;
    const arma::vec diag_Omega = arma::diagvec(Omega);
    const arma::mat ones_row   = arma::ones(n, 1);

    arma::mat S2   = S % S;
    arma::mat logS = arma::log(S);

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0;
    int last_status = 5;

    for (int it = 0; it < maxiter; it++) {
        S2 = S % S;
        arma::mat Z = O + X * B + M;
        arma::mat A = arma::exp(Z + 0.5 * S2);

        // ---- Diagonal Newton step for B ----
        arma::mat grad_B = Xw.t() * (A - Y);
        arma::mat hess_B = Xw2.t() * A;
        hess_B.clamp(1e-10, arma::datum::inf);
        arma::mat step_B  = grad_B / hess_B;
        arma::mat XstepB  = X * step_B;
        double f0_B   = arma::accu(w.t() * (A - Y % Z));
        double slope_B = -arma::accu(grad_B % step_B);
        double alpha_B = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            arma::mat Zt = Z - alpha_B * XstepB;
            if (arma::accu(w.t() * (arma::exp(Zt + 0.5*S2) - Y % Zt))
                <= f0_B + c1 * alpha_B * slope_B) break;
            alpha_B *= 0.5;
        }
        B -= alpha_B * step_B;
        Z  = O + X * B + M;
        A  = arma::exp(Z + 0.5 * S2);

        // ---- Diagonal Newton step for M ----
        arma::mat MO     = M * Omega;
        arma::mat grad_M = MO + A - Y; grad_M.each_col() %= w;
        arma::mat hess_M = A + ones_row * diag_Omega.t();
        hess_M.each_col() %= w;
        hess_M.clamp(1e-10, arma::datum::inf);
        arma::mat step_M = grad_M / hess_M;
        double f0_M   = arma::accu(w.t() * (A - Y % Z))
                      + 0.5 * arma::accu(MO % (M.each_col() % w));
        double slope_M = -arma::accu(grad_M % step_M);
        double alpha_M = 1.0;
        for (int ls = 0; ls < 20; ls++) {
            arma::mat Mt  = M - alpha_M * step_M;
            arma::mat Zt  = Z - alpha_M * step_M;
            arma::mat At  = arma::exp(Zt + 0.5 * S2);
            arma::mat MOt = Mt * Omega;
            if (arma::accu(w.t() * (At - Y % Zt))
                + 0.5 * arma::accu(MOt % (Mt.each_col() % w))
                <= f0_M + c1 * alpha_M * slope_M) break;
            alpha_M *= 0.5;
        }
        M  -= alpha_M * step_M;
        MO  = M * Omega;
        Z   = O + X * B + M;

        // ---- Fixed-point update for logS (overflow-safe) ----
        A = arma::exp(Z + 0.5 * S2);
        {
            const arma::mat logS_cand = -0.5 * arma::log(A + ones_row * diag_Omega.t());
            const arma::mat logS_ub   = 0.5 * arma::log(arma::clamp(700. - Z, 1., arma::datum::inf));
            logS = arma::clamp(arma::min(logS_cand, logS_ub), -20., arma::datum::inf);
        }
        S  = arma::exp(logS);
        S2 = S % S;

        // ---- Objective for convergence ----
        A = arma::exp(Z + 0.5 * S2);
        double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * arma::trunc_log(S2)))
                   + 0.5 * (arma::accu(MO % (M.each_col() % w))
                            + arma::dot(diag_Omega, (w.t() * S2).t()));
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && std::abs(obj - obj_prev) < ftol * (1.0 + std::abs(obj_prev))) {
            last_status = 3; break;
        }
        obj_prev = obj;
    }

    // ---- Final output ----
    S2 = S % S;
    arma::mat Sigma = (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2)) / w_bar;
    arma::mat Z = O + X * B + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec loglik = arma::sum(Y % Z - A - 0.5 * ((M * Omega) % M - arma::log(S2) + S2 * arma::diagmat(Omega)), 1)
                     + 0.5 * std::real(arma::log_det(Omega)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B",     B    ),
        Rcpp::Named("M",     M    ),
        Rcpp::Named("S",     S    ),
        Rcpp::Named("Z",     Z    ),
        Rcpp::Named("A",     A    ),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji",    Ji   ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status   ),
            Rcpp::Named("backend",    "newton"      ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", total_iter    )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// Fixed inverse covariance (Omega)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_fixed(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]);     // (d,p)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]);     // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(params["S"]);     // (n,p)
    const auto  Omega = Rcpp::as<arma::mat>(params["Omega"]); // covinv (p,p)

    const auto metadata = tuple_metadata(init_B, init_M, init_S);
    enum { B_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;

    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &Omega, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * trace(Omega * nSigma);

        metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);
        metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1)) ;

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat Sigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) / accu(w);
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * ((M * Omega) % M - log(S2) + S2 * diagmat(Omega)), 1) +
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
          Rcpp::Named("status", static_cast<int>(result.status)),
          Rcpp::Named("backend", "nlopt"),
          Rcpp::Named("objective", objective_vec),
          Rcpp::Named("iterations", result.nb_iterations)
      ))
    );
}

