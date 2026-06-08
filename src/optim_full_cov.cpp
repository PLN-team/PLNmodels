#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Fully parametrized covariance

// [[Rcpp::export]]
Rcpp::List nlopt_optimize(
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
            const arma::mat B = metadata.map<B_ID>(params);
            const arma::mat M = metadata.map<M_ID>(params);
            const arma::mat S = metadata.map<S_ID>(params);
            arma::mat S2 = S % S;
            arma::mat Z = O + X * B + M;
            arma::mat A = exp(Z + 0.5 * S2);
            arma::mat MO = M * Omega;           // cached: reused in objective and M-gradient
            const arma::rowvec wS2 = w.t() * S2;

            double objective = accu(w.t() * (A - Y % Z - 0.5 * trunc_log(S2)))
                             + 0.5 * (accu(MO % (M.each_col() % w)) + dot(Omega_diag, wS2.t()));

            metadata.map<B_ID>(grad) = Xw.t() * (A - Y);
            metadata.map<M_ID>(grad) = diagmat(w) * (MO + A - Y);
            metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % Omega_diag.t() + S % A - pow(S, -1));

            objective_vec.push_back(objective);
            return objective;
        };

        OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
        total_iterations += result.nb_iterations;
        last_status = static_cast<int>(result.status);

        // M-step: update Omega analytically (one inv_sympd per EM iteration)
        arma::mat M = metadata.copy<M_ID>(parameters.data());
        arma::mat S = metadata.copy<S_ID>(parameters.data());
        arma::mat S2 = S % S;
        arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
        Omega = inv_sympd(Sigma);

        // ELBO after M-step: trace(Omega*nSigma) = w_bar*p, log_det(Omega) = -log_det(Sigma)
        arma::mat B = metadata.copy<B_ID>(parameters.data());
        arma::mat Z = O + X * B + M;
        arma::mat A = exp(Z + 0.5 * S2);
        double elbo = accu(w.t() * (Y % Z - A + 0.5 * trunc_log(S2)))
                    - 0.5 * w_bar * real(log_det(Sigma));

        if (em_iter > 0 && converged(elbo, elbo_prev, em_ftol)) break;
        elbo_prev = elbo;
    }

    // Final extraction
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    // Omega already updated from the last M-step
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * log(S2) - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) +
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
// Full covariance PLN — coordinate-Newton optimizer

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_newton(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // (n,d)
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // (n,p)
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // (n)
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);        // (d,p)
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);        // (n,p)
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);        // (n,p), >0

    const int maxiter  = config.containsElementNamed("maxeval")     ? Rcpp::as<int>(config["maxeval"])        : 200;
    const double ftol  = config.containsElementNamed("ftol_rel")    ? Rcpp::as<double>(config["ftol_rel"])    : 1e-8;
    const int max_em   = config.containsElementNamed("max_em_iter") ? Rcpp::as<int>(config["max_em_iter"])    : 50;
    const double em_tol= config.containsElementNamed("em_ftol")     ? Rcpp::as<double>(config["em_ftol"])     : 1e-8;

    const int n = Y.n_rows;
    const double w_bar = arma::accu(w);
    const double c1    = 1e-4;

    // Precompute weighted X matrices (fixed throughout)
    const arma::mat Xw  = X.each_col() % w;            // (n,d)  X_ik * w_i
    arma::mat Xw2 = X % X; Xw2.each_col() %= w;        // (n,d)  X_ik^2 * w_i

    // Initial Omega from starting M, S
    arma::mat S2    = S % S;
    arma::mat Sigma = (1./w_bar) * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2));
    arma::mat Omega = arma::inv_sympd(Sigma);
    arma::mat logS  = arma::log(S);

    std::vector<double> objective_vec;
    double elbo_prev  = -arma::datum::inf;
    int    total_iter = 0;
    int    last_status = 5; // maxiter reached by default

    for (int em = 0; em < max_em; em++) {
        const arma::vec diag_Omega = arma::diagvec(Omega);
        const arma::mat ones_row   = arma::ones(n, 1);

        double obj_prev = arma::datum::inf;

        for (int it = 0; it < maxiter; it++) {
            S2 = S % S;
            arma::mat Z = O + X * B + M;
            arma::mat A = arma::exp(Z + 0.5 * S2);

            // ---- Diagonal Newton step for B ----
            newton_step_B(Xw, Xw2, X, Y, O, w, M, S2, B, Z, A);

            // ---- Diagonal Newton step for M ----
            arma::mat MO     = M * Omega;                              // (n,p)
            arma::mat grad_M = MO + A - Y; grad_M.each_col() %= w;    // (n,p)
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
            fixed_point_logS(logS, S, S2, Z, A, ones_row * diag_Omega.t());

            // ---- Objective for inner convergence ----
            A = arma::exp(Z + 0.5 * S2);
            double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * arma::trunc_log(S2)))
                       + 0.5 * (arma::accu(MO % (M.each_col() % w))
                                + arma::dot(diag_Omega, (w.t() * S2).t()));
            objective_vec.push_back(obj);
            total_iter++;

            if (it > 0 && converged(obj, obj_prev, ftol)) { last_status = 3; break; }
            obj_prev = obj;
        }

        // ---- M-step: update Omega analytically ----
        S2    = S % S;
        Sigma = (1./w_bar) * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2));
        Omega = arma::inv_sympd(Sigma);

        // ---- Outer EM convergence on ELBO ----
        arma::mat Z    = O + X * B + M;
        arma::mat A    = arma::exp(Z + 0.5 * S2);
        double elbo = arma::accu(w.t() * (Y % Z - A + 0.5 * arma::trunc_log(S2)))
                    - 0.5 * w_bar * real(arma::log_det(Sigma));
        if (em > 0 && converged(elbo, elbo_prev, em_tol)) { last_status = 3; break; }
        elbo_prev = elbo;
    }

    // ---- Final output ----
    S2 = S % S;
    arma::mat Z = O + X * B + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec loglik = arma::sum(Y % Z - A + 0.5 * arma::log(S2)
                         - 0.5 * ((M * Omega) % M + S2 * arma::diagmat(Omega)), 1)
                     + 0.5 * real(arma::log_det(Omega)) + ki(Y);

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
// VE full — coordinate-Newton (M and S only, B and Omega fixed)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_newton(
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

    const int maxiter = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 200;
    const double ftol = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-8;

    const int n = Y.n_rows;
    const double c1 = 1e-4;
    const arma::vec diag_Omega = arma::diagvec(Omega);
    const arma::mat ones_row   = arma::ones(n, 1);
    const arma::mat XB         = X * B;

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
        Z   = O + XB + M;

        // ---- Fixed-point update for logS ----
        A    = arma::exp(Z + 0.5 * S2);
        logS = arma::clamp(-0.5 * arma::log(A + ones_row * diag_Omega.t()), -20., 10.);
        S    = arma::exp(logS);
        S2   = S % S;

        // ---- Objective for convergence check ----
        A = arma::exp(Z + 0.5 * S2);
        double obj = arma::accu(w.t() * (A - Y % Z - 0.5 * arma::trunc_log(S2)))
                   + 0.5 * (arma::accu(MO % (M.each_col() % w))
                            + arma::dot(diag_Omega, (w.t() * S2).t()));
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
    }

    // ---- Final output ----
    S2 = S % S;
    arma::mat Z = O + XB + M;
    arma::mat A = arma::exp(Z + 0.5 * S2);
    arma::vec loglik = arma::sum(Y % Z - A + 0.5 * arma::log(S2)
                         - 0.5 * ((M * Omega) % M + S2 * arma::diagmat(Omega)), 1)
                     + 0.5 * real(arma::log_det(Omega)) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S")  = S,
        Rcpp::Named("Ji") = Ji,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     3          ),
            Rcpp::Named("backend",    "newton"   ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE full

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep(
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
    metadata.map<S_ID>(parameters.data()) = init_S;

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec ;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat XB_vestep   = X * B;
    const arma::vec Omega_diag_v = diagvec(Omega);

    auto objective_and_grad = [&metadata, &O, &XB_vestep, &Y, &w, &Omega, &Omega_diag_v, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + XB_vestep + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat MO = M * Omega;
        const arma::rowvec wS2 = w.t() * S2;
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2)))
                         + 0.5 * (accu(MO % (M.each_col() % w)) + dot(Omega_diag_v, wS2.t()));

        metadata.map<M_ID>(grad) = diagmat(w) * (MO + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % Omega_diag_v.t() + S % A - pow(S, -1));

        objective_vec.push_back(objective) ;

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Element-wise log-likelihood
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * log(S2) - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) +
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
