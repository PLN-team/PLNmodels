#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

// ---------------------------------------------------------------------------------------
// Rank-constrained PLN — Spectral Gradient Method (Barzilai-Borwein + GLL/Armijo)
//
// Parameters: B (d,p), C (p,q), M (n,q), S exact via ψ = −log(1 + A·C²)
// Z = O + X·B + M·C',   A = exp(Z + ½·S²·C²ᵀ)
//
// Per-element BB step sizes (equivalent to NLOPT CCSAQ's σ):
//   σᵢ = |Δgᵢ| / |sᵢ|   for each element i of B, C, M
//
// Line search: GLL nonmonotone Armijo — accepts if
//   f_new ≤ max_{j=0..M-1} f_{k-j} + c₁·scale·slope
// GLL reduces backtracking frequency vs monotone Armijo (the nonmonotone window
// allows temporary increases, so the BB step is accepted more often on the first try).
//
// Convergence: consecutive ftol (default 1e-9 in config_default_spectral).
// ftol=1e-8 is too loose — GLL oscillations are O(1e-8·|f|) and fool the
// criterion into stopping prematurely at suboptimal minima.  1e-9 gives
// quality matching or exceeding NLOPT CCSAQ at ≤2× wall-clock time.

static inline double safe_ratio(double a, double b, double fallback) {
    return (b > 1e-20) ? std::clamp(a / b, 1e-12, 1e12) : fallback;
}

// Per-element σ update: σ = |y| / |s|, fall back to current value when |s| is tiny
static void update_sigma(arma::mat & sigma,
                         const arma::mat & y_abs,
                         const arma::mat & s_abs)
{
    for (arma::uword i = 0; i < sigma.n_elem; i++)
        sigma(i) = safe_ratio(y_abs(i), s_abs(i), sigma(i));
}

// [[Rcpp::export]]
Rcpp::List spectral_optimize_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B = Rcpp::as<arma::mat>(params["B"]);
    arma::mat C = Rcpp::as<arma::mat>(params["C"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter   = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 10000;
    const double ftol      = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-10;
    const double c_armijo  = 1e-4;
    const int    nm_window = 10;   // GLL nonmonotone window size

    const arma::mat Xw = X.each_col() % w;

    // Initialise geometry with S exact
    arma::mat C2  = C % C;
    arma::mat psi = arma::log(S % S);
    arma::mat S2  = arma::exp(psi);
    arma::mat Z   = O + X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    // Per-element BB σ (same layout as B, C, M)
    arma::mat sigma_B(arma::size(B), arma::fill::ones);
    arma::mat sigma_C(arma::size(C), arma::fill::ones);
    arma::mat sigma_M(arma::size(M), arma::fill::ones);
    arma::mat sB_prev, sC_prev, sM_prev;
    arma::mat gB_prev, gC_prev, gM_prev;
    bool prev_valid = false;

    // GLL nonmonotone buffer
    std::deque<double> obj_buffer;

    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0, last_status = 5;
    const int win = 100;  // window for GLL-oscillation-robust convergence check

    for (int it = 0; it < maxiter; it++) {

        // ---- Objective ----
        double obj = arma::accu(w.t() * (A - Y % Z))
                   + 0.5 * arma::accu(w.t() * (M % M + S2 - psi - 1.));
        objective_vec.push_back(obj);
        total_iter++;
        // Primary: consecutive convergence (cheap, exact for monotone steps)
        if (it > 0 && converged(obj, obj_prev, ftol)) { last_status = 3; break; }
        obj_prev = obj;
        // Secondary: window-minimum stagnation every `win` iterations.
        // Compares min of last `win` evals to min of the preceding `win` evals.
        // Robust to GLL oscillations, which can delay consecutive convergence.
        if (it > 0 && it % win == 0 && (int)objective_vec.size() >= 2*win) {
            double m1 = *std::min_element(objective_vec.end()-win, objective_vec.end());
            double m2 = *std::min_element(objective_vec.end()-2*win, objective_vec.end()-win);
            if (converged(m1, m2, ftol)) { last_status = 3; break; }
        }

        // Update GLL window and reference
        obj_buffer.push_back(obj);
        if ((int)obj_buffer.size() > nm_window) obj_buffer.pop_front();
        double f_ref = *std::max_element(obj_buffer.begin(), obj_buffer.end());

        // ---- Joint gradient ∇(B, C, M) ----
        arma::mat AmY = A - Y;  AmY.each_col() %= w;
        arma::mat gB  = Xw.t() * (A - Y);
        arma::mat gC  = AmY.t() * M + (A.t() * (S2.each_col() % w)) % C;
        arma::mat gM  = (A - Y) * C + M;  gM.each_col() %= w;

        // ---- Per-element σ update (or initialise on first iter) ----
        if (prev_valid) {
            update_sigma(sigma_B, arma::abs(gB - gB_prev), arma::abs(sB_prev));
            update_sigma(sigma_C, arma::abs(gC - gC_prev), arma::abs(sC_prev));
            update_sigma(sigma_M, arma::abs(gM - gM_prev), arma::abs(sM_prev));
        } else {
            sigma_B = arma::clamp(arma::abs(gB), 1e-4, 1e12);
            sigma_C = arma::clamp(arma::abs(gC), 1e-4, 1e12);
            sigma_M = arma::clamp(arma::abs(gM), 1e-4, 1e12);
        }

        arma::mat dB = gB / sigma_B;
        arma::mat dC = gC / sigma_C;
        arma::mat dM = gM / sigma_M;
        double slope = -(arma::accu(dB % gB) + arma::accu(dC % gC) + arma::accu(dM % gM));

        // ---- GLL nonmonotone Armijo with S exact update ----
        double scale = 1.0;
        arma::mat B_t, C_t, M_t, C2_t, Z_t, A_t, psi_t, S2_t;
        double obj_t = arma::datum::inf;

        for (int ls = 0; ls < 20; ls++) {
            B_t  = B  - scale * dB;
            C_t  = C  - scale * dC;
            M_t  = M  - scale * dM;
            C2_t = C_t % C_t;
            Z_t  = O + X * B_t + M_t * C_t.t();
            A_t   = arma::exp(Z_t + 0.5 * S2 * C2_t.t());
            psi_t = arma::clamp(-arma::log(1. + A_t * C2_t), -40., 0.);
            S2_t  = arma::exp(psi_t);
            A_t   = arma::exp(Z_t + 0.5 * S2_t * C2_t.t());
            obj_t = arma::accu(w.t() * (A_t - Y % Z_t))
                  + 0.5 * arma::accu(w.t() * (M_t % M_t + S2_t - psi_t - 1.));
            if (obj_t <= f_ref + c_armijo * scale * slope) break;
            scale *= 0.5;
        }

        // σ scale-up when Armijo backtracked heavily — prevents the next step
        // from overshooting in the same direction
        if (scale < 0.125) {
            sigma_B /= scale;
            sigma_C /= scale;
            sigma_M /= scale;
        }

        sB_prev = B_t - B;  sC_prev = C_t - C;  sM_prev = M_t - M;
        gB_prev = gB;  gC_prev = gC;  gM_prev = gM;
        prev_valid = true;

        B = B_t;  C = C_t;  M = M_t;
        C2 = C2_t;  Z = Z_t;  A = A_t;  psi = psi_t;  S2 = S2_t;
    }

    // ---- Final output ----
    S = arma::exp(0.5 * psi);
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
            Rcpp::Named("backend",    "spectral"   ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE step: B and C fixed, update M (per-element BB + GLL) and S (exact)

// [[Rcpp::export]]
Rcpp::List spectral_optimize_vestep_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const arma::mat & B,
    const arma::mat & C,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat M = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S = Rcpp::as<arma::mat>(params["S"]);

    const int    maxiter   = config.containsElementNamed("maxeval")  ? Rcpp::as<int>(config["maxeval"])     : 10000;
    const double ftol      = config.containsElementNamed("ftol_rel") ? Rcpp::as<double>(config["ftol_rel"]) : 1e-10;
    const double c_armijo  = 1e-4;
    const int    nm_window = 10;

    const arma::mat C2 = C % C;
    const arma::mat XB = X * B;

    arma::mat psi = arma::log(S % S);
    arma::mat S2  = arma::exp(psi);
    arma::mat Z   = O + XB + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 0.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    arma::mat sigma_M(arma::size(M), arma::fill::ones);
    arma::mat sM_prev, gM_prev;
    bool prev_valid = false;

    std::deque<double> obj_buffer;
    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0;
    const int win = 100;

    for (int it = 0; it < maxiter; it++) {
        double obj = arma::accu(w.t() * (A - Y % Z))
                   + 0.5 * arma::accu(w.t() * (M % M + S2 - psi - 1.));
        objective_vec.push_back(obj);
        total_iter++;
        if (it > 0 && converged(obj, obj_prev, ftol)) break;
        obj_prev = obj;
        if (it > 0 && it % win == 0 && (int)objective_vec.size() >= 2*win) {
            double m1 = *std::min_element(objective_vec.end()-win, objective_vec.end());
            double m2 = *std::min_element(objective_vec.end()-2*win, objective_vec.end()-win);
            if (converged(m1, m2, ftol)) break;
        }

        obj_buffer.push_back(obj);
        if ((int)obj_buffer.size() > nm_window) obj_buffer.pop_front();
        double f_ref = *std::max_element(obj_buffer.begin(), obj_buffer.end());

        arma::mat gM = (A - Y) * C + M;  gM.each_col() %= w;

        if (prev_valid) {
            update_sigma(sigma_M, arma::abs(gM - gM_prev), arma::abs(sM_prev));
        } else {
            sigma_M = arma::clamp(arma::abs(gM), 1e-4, 1e12);
        }

        arma::mat dM = gM / sigma_M;
        double slope = -arma::accu(dM % gM);
        double scale = 1.0;
        arma::mat M_t, Z_t, A_t, psi_t, S2_t;
        double obj_t = arma::datum::inf;

        for (int ls = 0; ls < 20; ls++) {
            M_t   = M - scale * dM;
            Z_t   = O + XB + M_t * C.t();
            A_t   = arma::exp(Z_t + 0.5 * S2 * C2.t());
            psi_t = arma::clamp(-arma::log(1. + A_t * C2), -40., 0.);
            S2_t  = arma::exp(psi_t);
            A_t   = arma::exp(Z_t + 0.5 * S2_t * C2.t());
            obj_t = arma::accu(w.t() * (A_t - Y % Z_t))
                  + 0.5 * arma::accu(w.t() * (M_t % M_t + S2_t - psi_t - 1.));
            if (obj_t <= f_ref + c_armijo * scale * slope) break;
            scale *= 0.5;
        }

        if (scale < 0.125) sigma_M /= scale;

        sM_prev = M_t - M;
        gM_prev = gM;
        prev_valid = true;

        M = M_t;  Z = Z_t;  A = A_t;  psi = psi_t;  S2 = S2_t;
    }

    S = arma::exp(0.5 * psi);
    Z  = O + XB + M * C.t();
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
            Rcpp::Named("status",     3            ),
            Rcpp::Named("backend",    "spectral"   ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}
