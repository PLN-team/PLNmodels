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

    const arma::uword n = Y.n_rows;
    const arma::uword q = C.n_cols;

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

        // Precompute once per iteration — reused in B and M Armijo (saves O(n*q*p) per backtrack)
        arma::mat half_S2C2t = 0.5 * S2 * C2.t();

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
                if (arma::accu(w.t() * (arma::exp(Zt + half_S2C2t) - Y % Zt))
                    <= f0_B + c1 * alpha_B * slope_B) break;
                alpha_B *= 0.5;
            }
            B -= alpha_B * step_B;
            Z  = O + X * B + M * C.t();
            A  = arma::exp(Z + half_S2C2t);
        }

        // ---- M step (block q×q Newton per observation + Armijo) ----
        // Exact Hessian per row i: Hᵢ = Cᵀ diag(Aᵢ) C + Iq  (q×q, always PD)
        {
            arma::mat grad_M = (A - Y) * C + M;  // n×q, unweighted gradient
            arma::mat step_M(n, q);
            for (int i = 0; i < n; i++) {
                // Scale rows of C by A_{ij}: CAi_{jk} = C_{jk} * A_{ij}
                arma::mat CAi = C.each_col() % A.row(i).t();  // p×q, O(p*q)
                arma::mat Hi  = CAi.t() * C;                   // q×q = C'diag(Aᵢ)C, O(p*q²)
                Hi.diag() += 1.0;
                step_M.row(i) = arma::solve(Hi, grad_M.row(i).t(), arma::solve_opts::fast).t();
            }
            arma::mat stepMCt  = step_M * C.t();
            arma::mat grad_M_w = grad_M;  grad_M_w.each_col() %= w;  // for slope
            double f0_M    = arma::accu(w.t() * (A - Y % Z))
                           + 0.5 * arma::as_scalar(w.t() * arma::sum(M % M, 1));
            double slope_M = -arma::accu(grad_M_w % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt = M - alpha_M * step_M;
                arma::mat Zt = Z - alpha_M * stepMCt;
                arma::mat At = arma::exp(Zt + half_S2C2t);
                if (arma::accu(w.t() * (At - Y % Zt))
                    + 0.5 * arma::as_scalar(w.t() * arma::sum(Mt % Mt, 1))
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + X * B + M * C.t();
            A  = arma::exp(Z + half_S2C2t);
        }

        // ---- C step (block q×q Newton per gene + Armijo) ----
        // Exact Hessian per gene j: Hⱼ = F̃ⱼᵀ diag(w⊙Aⱼ) F̃ⱼ + diag(S2ᵀ(w⊙Aⱼ))
        // where F̃ⱼ_{ik} = Mᵢₖ + S²ᵢₖ Cⱼₖ  (avoids n×p temporaries of diagonal version)
        {
            const arma::uword p = C.n_rows;
            const arma::mat WA   = A.each_col() % w;          // n×p
            arma::mat AmY        = A - Y;  AmY.each_col() %= w;
            arma::mat grad_C     = AmY.t() * M + (WA.t() * S2) % C;  // p×q
            arma::mat step_C(p, q);
            for (arma::uword j = 0; j < p; j++) {
                // F̃ⱼ = M + S2 .* C[j,:]  (broadcast C.row(j) over n rows)
                arma::mat Fj    = M + S2.each_row() % C.row(j);       // n×q, O(n*q)
                arma::mat Fj_sc = Fj.each_col() % WA.col(j);          // n×q, scale by wᵢAᵢⱼ
                arma::mat Hj    = Fj_sc.t() * Fj;                     // q×q, O(n*q²)
                Hj.diag()      += S2.t() * WA.col(j);                 // + diag(S2ᵀ wA_j)
                step_C.row(j)   = arma::solve(Hj, grad_C.row(j).t(), arma::solve_opts::fast).t();
            }
            const arma::mat Mstep_Ct = M * step_C.t();
            const arma::mat S2_CsCt  = S2 * (C % step_C).t();
            const arma::mat S2_sC2t  = S2 * (step_C % step_C).t();
            double f0_C    = arma::accu(w.t() * (A - Y % Z));
            double slope_C = -arma::accu(grad_C % step_C);
            double alpha_C = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Zt    = Z - alpha_C * Mstep_Ct;
                arma::mat half_t = half_S2C2t - alpha_C * S2_CsCt
                                 + (0.5 * alpha_C * alpha_C) * S2_sC2t;
                arma::mat At    = arma::exp(Zt + half_t);
                if (arma::accu(w.t() * (At - Y % Zt))
                    <= f0_C + c1 * alpha_C * slope_C) break;
                alpha_C *= 0.5;
            }
            C  -= alpha_C * step_C;
            C2  = C % C;
            Z   = O + X * B + M * C.t();
            half_S2C2t = 0.5 * S2 * C2.t();
            A   = arma::exp(Z + half_S2C2t);
        }

        // ---- S step (exact minimiser in ψ = logS² space) ----
        // ∂F/∂ψᵢₖ = 0  ⟹  ψᵢₖ = −log(1 + (A·C²)ᵢₖ)
        // KL term: S² − logS² − 1 = exp(ψ) − ψ − 1
        {
            psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
            S2  = arma::exp(psi);
            A   = arma::exp(Z + 0.5 * S2 * C2.t());   // S2 changed — full recompute
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

    const arma::uword n = Y.n_rows;
    const arma::uword q = C.n_cols;

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

        const arma::mat half_S2C2t = 0.5 * S2 * C2.t();

        // ---- M step (block q×q Newton per observation + Armijo) ----
        {
            arma::mat grad_M = (A - Y) * C + M;  // n×q, unweighted
            arma::mat step_M(n, q);
            for (int i = 0; i < n; i++) {
                arma::mat CAi = C.each_col() % A.row(i).t();
                arma::mat Hi  = CAi.t() * C;
                Hi.diag() += 1.0;
                step_M.row(i) = arma::solve(Hi, grad_M.row(i).t(), arma::solve_opts::fast).t();
            }
            arma::mat stepMCt  = step_M * C.t();
            arma::mat grad_M_w = grad_M;  grad_M_w.each_col() %= w;
            double f0_M    = arma::accu(w.t() * (A - Y % Z))
                           + 0.5 * arma::as_scalar(w.t() * arma::sum(M % M, 1));
            double slope_M = -arma::accu(grad_M_w % step_M);
            double alpha_M = 1.0;
            for (int ls = 0; ls < 20; ls++) {
                arma::mat Mt = M - alpha_M * step_M;
                arma::mat Zt = Z - alpha_M * stepMCt;
                arma::mat At = arma::exp(Zt + half_S2C2t);
                if (arma::accu(w.t() * (At - Y % Zt))
                    + 0.5 * arma::as_scalar(w.t() * arma::sum(Mt % Mt, 1))
                    <= f0_M + c1 * alpha_M * slope_M) break;
                alpha_M *= 0.5;
            }
            M -= alpha_M * step_M;
            Z  = O + XB + M * C.t();
            A  = arma::exp(Z + half_S2C2t);
        }

        // ---- S step (exact minimiser in ψ = logS² space) ----
        {
            psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
            S2  = arma::exp(psi);
            A   = arma::exp(Z + 0.5 * S2 * C2.t());   // S2 changed — full recompute
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
