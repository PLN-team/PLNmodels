#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

// ---------------------------------------------------------------------------------------
// Rank-constrained covariance PLN — coordinate-Newton optimizer
//
// Parameters: B (d,p), C (p,q), M (n,q), S (n,q)
// Z = O + X*B + M*C',   A = exp(Z + 0.5 * S²*C²ᵀ)
// Objective: Σᵢ wᵢ (A - Y⊙Z) + ½ Σᵢ wᵢ Σₖ (Mᵢₖ² + Sᵢₖ² − log Sᵢₖ² − 1)
//
// Update order per iteration: B → C → M → S
//   B, M : diagonal Newton + Armijo (standard)
//   C    : diagonal Newton (Gauss-Newton Hessian) + Armijo (full Z and A recomputed)
//   S    : closed-form exact minimiser in ψ = logS² space: ψ = −log(1 + A·C²)

// [[Rcpp::export]]
Rcpp::List newton_optimize_rank(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, C, M, S)
    const Rcpp::List & config  // List of config values
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);
    arma::mat C = Rcpp::as<arma::mat>(params["C"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 200;
    const double ftol    = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-8;

    const double c1 = 1e-4;
    const arma::mat Xw  = X.each_col() % w;
    arma::mat Xw2 = X % X; Xw2.each_col() %= w;

    arma::mat C2  = C % C;
    arma::mat psi = arma::log(S % S);        // ψ = log(S²), work in logS² space throughout
    arma::mat S2  = arma::exp(psi);
    arma::mat Z   = O + X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter  = 0;
    int last_status = 5;

    for (int it = 0; it < maxiter; it++) {

        // ---- B step (diagonal Newton + Armijo) ----
        {
            arma::mat grad_B  = Xw.t() * (A - Y);
            arma::mat hess_B  = Xw2.t() * A;
            hess_B.clamp(1e-10, arma::datum::inf);
            arma::mat step_B  = grad_B / hess_B;
            arma::mat XstepB  = X * step_B;
            double    f0_B    = arma::accu(w.t() * (A - Y % Z));
            double    slope_B = -arma::accu(grad_B % step_B);
            double    alpha_B = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Zt = Z - alpha_B * XstepB;
                if (arma::accu(w.t() * (arma::exp(Zt + 0.5 * S2 * C2.t()) - Y % Zt))
                    <= f0_B + c1 * alpha_B * slope_B) break;
                alpha_B *= 0.5;
            }
            B -= alpha_B * step_B;
            Z  = O + X * B + M * C.t();
            A  = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- M step (diagonal Newton + Armijo) ----
        // grad_M = diag(w) * ((A-Y)*C + M),  hess_M = diag(w) * (A*C² + 1)
        {
            arma::mat grad_M  = (A - Y) * C + M;  grad_M.each_col() %= w;
            arma::mat hess_M  = A * C2 + 1.;      hess_M.each_col() %= w;
            hess_M.clamp(1e-10, arma::datum::inf);
            arma::mat step_M  = grad_M / hess_M;
            arma::mat stepMCt = step_M * C.t();
            double    f0_M    = arma::accu(w.t() * (A - Y % Z))
                              + 0.5 * arma::as_scalar(w.t() * arma::sum(M % M, 1));
            double    slope_M = -arma::accu(grad_M % step_M);
            double    alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt = M - alpha_M * step_M;
                arma::mat Zt = Z - alpha_M * stepMCt;
                arma::mat At = arma::exp(Zt + 0.5 * S2 * C2.t());
                if (arma::accu(w.t() * (At - Y % Zt))
                    + 0.5 * arma::as_scalar(w.t() * arma::sum(Mt % Mt, 1))
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + X * B + M * C.t();
            A  = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- C step (exact diagonal Newton + Armijo) ----
        // Exact diagonal Hessian: h_{jk} = Σᵢ wᵢ Aᵢⱼ [(Mᵢₖ + S²ᵢₖ Cⱼₖ)² + S²ᵢₖ]
        {
            const arma::uword p = C.n_rows, q = C.n_cols;
            arma::mat AmY = A - Y;  AmY.each_col() %= w;
            arma::mat grad_C  = AmY.t() * M + (A.t() * (S2.each_col() % w)) % C;
            arma::mat hess_C(p, q);
            for (arma::uword k = 0; k < q; k++) {
                // Fk_{ij} = M_{ik} + S2_{ik} * C_{jk}  (n×p outer product)
                arma::mat Fk = M.col(k) * arma::ones(1, p)
                             + S2.col(k) * C.col(k).t();
                hess_C.col(k) = (A % (arma::square(Fk) + S2.col(k) * arma::ones(1, p))).t() * w;
            }
            hess_C.clamp(1e-10, arma::datum::inf);
            arma::mat step_C  = grad_C / hess_C;
            double    f0_C    = arma::accu(w.t() * (A - Y % Z));
            double    slope_C = -arma::accu(grad_C % step_C);
            double    alpha_C = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Ct = C  - alpha_C * step_C;
                arma::mat Zt = O + X * B + M * Ct.t();
                arma::mat At = arma::exp(Zt + 0.5 * S2 * (Ct % Ct).t());
                if (arma::accu(w.t() * (At - Y % Zt))
                    <= f0_C + c1 * alpha_C * slope_C) break;
                alpha_C *= 0.5;
            }
            C  -= alpha_C * step_C;
            C2  = C % C;
            Z   = O + X * B + M * C.t();
            A   = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- S step (exact minimiser in ψ = logS² space) ----
        // ∂F/∂ψᵢₖ = 0  ⟹  ψᵢₖ = −log(1 + (A·C²)ᵢₖ)
        // KL term: S² − logS² − 1 = exp(ψ) − ψ − 1
        {
            psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
            S2  = arma::exp(psi);
            A   = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- Convergence check ----
        double obj = arma::accu(w.t() * (A - Y % Z))
                   + 0.5 * arma::accu(w.t() * (M % M + S2 - psi - 1.));
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) { last_status = 3; break; }
        obj_prev = obj;
    }

    // ---- Final output ----
    S2 = arma::exp(psi);
    S  = arma::exp(0.5 * psi);
    Z  = O + X * B + M * C.t();
    A  = arma::exp(Z + 0.5 * S2 * C2.t());
    const double w_bar = arma::accu(w);
    arma::mat nSig = M.t() * (M.each_col() % w) + arma::diagmat(arma::sum(S2.each_col() % w, 0));
    arma::mat Sigma = C * nSig * C.t() / w_bar;
    arma::mat Omega = C * arma::inv_sympd(nSig / w_bar) * C.t();
    arma::vec loglik = arma::sum(Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1)
                     + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B",     B    ),
        Rcpp::Named("C",     C    ),
        Rcpp::Named("M",     M    ),
        Rcpp::Named("S",     S    ),
        Rcpp::Named("Z",     Z    ),
        Rcpp::Named("A",     A    ),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji",    Ji   ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status  ),
            Rcpp::Named("backend",    "newton"     ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE rank — coordinate-Newton (M and S only, B and C fixed)

// [[Rcpp::export]]
Rcpp::List newton_optimize_vestep_rank(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S)
    const arma::mat & B,       // (d,p) fixed
    const arma::mat & C,       // (p,q) fixed
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

    const double  c1  = 1e-4;
    const arma::mat C2  = C % C;
    const arma::mat XB  = X * B;

    arma::mat psi = arma::log(S % S);   // ψ = log(S²)
    arma::mat S2  = arma::exp(psi);
    arma::mat Z   = O + XB + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int    total_iter = 0;

    for (int it = 0; it < maxiter; it++) {

        // ---- M step (diagonal Newton + Armijo) ----
        {
            arma::mat grad_M  = (A - Y) * C + M;  grad_M.each_col() %= w;
            arma::mat hess_M  = A * C2 + 1.;      hess_M.each_col() %= w;
            hess_M.clamp(1e-10, arma::datum::inf);
            arma::mat step_M  = grad_M / hess_M;
            arma::mat stepMCt = step_M * C.t();
            double    f0_M    = arma::accu(w.t() * (A - Y % Z))
                              + 0.5 * arma::as_scalar(w.t() * arma::sum(M % M, 1));
            double    slope_M = -arma::accu(grad_M % step_M);
            double    alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt = M - alpha_M * step_M;
                arma::mat Zt = Z - alpha_M * stepMCt;
                arma::mat At = arma::exp(Zt + 0.5 * S2 * C2.t());
                if (arma::accu(w.t() * (At - Y % Zt))
                    + 0.5 * arma::as_scalar(w.t() * arma::sum(Mt % Mt, 1))
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + XB + M * C.t();
            A  = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- S step (exact minimiser in ψ = logS² space) ----
        {
            psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
            S2  = arma::exp(psi);
            A   = arma::exp(Z + 0.5 * S2 * C2.t());
        }

        // ---- Convergence check ----
        double obj = arma::accu(w.t() * (A - Y % Z))
                   + 0.5 * arma::accu(w.t() * (M % M + S2 - psi - 1.));
        objective_vec.push_back(obj);
        total_iter++;

        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
    }

    // ---- Final output ----
    S2 = arma::exp(psi);
    S  = arma::exp(0.5 * psi);
    Z  = O + XB + M * C.t();
    A  = arma::exp(Z + 0.5 * S2 * C2.t());
    arma::vec loglik = arma::sum(Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1)
                     + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S")  = S,
        Rcpp::Named("Ji") = Ji,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     3           ),
            Rcpp::Named("backend",    "newton"    ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter  )
        ))
    );
}
