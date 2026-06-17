#pragma once

#include <RcppArmadillo.h>

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::mat logfact_mat(arma::mat y) {
    y.replace(0., 1.);
    return y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2.;
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
    return -logfact(std::move(y)) + 0.5 * double(p) ;
}

inline arma::mat logistic(arma::mat M) {
  arma::mat e = arma::trunc_exp(M);
  return e / (1. + e);
}

inline arma::mat logit(arma::mat M) {
  return arma::trunc_log(M) - arma::trunc_log(1 - M) ;
}

// ---- Relative convergence test: |val - prev| < tol * (1 + |prev|) ----
inline bool converged(double val, double prev, double tol) {
    return std::abs(val - prev) < tol * (1.0 + std::abs(prev));
}

// ---- Config extraction for builtin Newton optimizers ----
// Centralises the containsElementNamed pattern used by builtin_optim_pln.*, nlopt_optim_pln.h
// and builtin_optim_plnpca.cpp. Defaults are constructor arguments (not in-class initializers)
// so callers with different historical defaults (e.g. PLNPCA's 10000/1e-9) can opt in without
// changing behavior for the rest of the codebase.
struct NewtonConfig {
    int    maxiter;
    double ftol;
    int    max_em;
    double em_tol;
    explicit NewtonConfig(const Rcpp::List & cfg,
                           int default_maxiter = 200, double default_ftol = 1e-8,
                           int default_max_em = 50,    double default_em_tol = 1e-8)
        : maxiter(default_maxiter), ftol(default_ftol), max_em(default_max_em), em_tol(default_em_tol) {
        if (cfg.containsElementNamed("maxeval"))  maxiter = Rcpp::as<int>(cfg["maxeval"]);
        if (cfg.containsElementNamed("ftol_in"))  ftol    = Rcpp::as<double>(cfg["ftol_in"]);
        if (cfg.containsElementNamed("maxit_em")) max_em  = Rcpp::as<int>(cfg["maxit_em"]);
        if (cfg.containsElementNamed("ftol_em"))  em_tol  = Rcpp::as<double>(cfg["ftol_em"]);
    }
};

// ---- PLN data extraction: Y, X, O, w from the R-side `data` list ----
// Centralises the Rcpp::as<...> pattern replicated across all PLN optimizers (builtin and nlopt).
struct PlnData {
    arma::mat Y, X, O;
    arma::vec w;
    explicit PlnData(const Rcpp::List & data) :
        Y(Rcpp::as<arma::mat>(data["Y"])),
        X(Rcpp::as<arma::mat>(data["X"])),
        O(Rcpp::as<arma::mat>(data["O"])),
        w(Rcpp::as<arma::vec>(data["w"])) {}
};

// ---- Result assembly for PLN optimizers ----
// Centralises the two Rcpp::List::create shapes replicated across all PLN optimizers.
// "Full" result: joint (B, Sigma/Omega) optimizers (nlopt EM/profiled, builtin Newton).
inline Rcpp::List make_pln_result(
    const arma::mat & B, const arma::mat & M, const arma::mat & S2,
    const arma::mat & Z, const arma::mat & A,
    const Rcpp::List & cov_out,   // List(Sigma, Omega)
    const arma::vec & loglik,     // per-observation log-likelihood Ji (auto-wrapped by Rcpp::Named)
    int status, const char * backend,
    const std::vector<double> & objective_vec, int iterations
) {
    return Rcpp::List::create(
        Rcpp::Named("B",     B               ),
        Rcpp::Named("M",     M               ),
        Rcpp::Named("S2",    S2              ),
        Rcpp::Named("Z",     Z               ),
        Rcpp::Named("A",     A               ),
        Rcpp::Named("Sigma", cov_out["Sigma"]),
        Rcpp::Named("Omega", cov_out["Omega"]),
        Rcpp::Named("Ji",    loglik          ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     status        ),
            Rcpp::Named("backend",    backend       ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", iterations    )
        ))
    );
}

// "Rank" result: like the "Full" shape above, plus the loadings matrix C (PLNPCA only).
inline Rcpp::List make_plnpca_result(
    const arma::mat & B, const arma::mat & C, const arma::mat & M, const arma::mat & S2,
    const arma::mat & Z, const arma::mat & A,
    const Rcpp::List & cov_out,   // List(Sigma, Omega)
    const arma::vec & loglik,
    int status, const char * backend,
    const std::vector<double> & objective_vec, int iterations
) {
    return Rcpp::List::create(
        Rcpp::Named("B",     B               ),
        Rcpp::Named("C",     C               ),
        Rcpp::Named("M",     M               ),
        Rcpp::Named("S2",    S2              ),
        Rcpp::Named("Z",     Z               ),
        Rcpp::Named("A",     A               ),
        Rcpp::Named("Sigma", cov_out["Sigma"]),
        Rcpp::Named("Omega", cov_out["Omega"]),
        Rcpp::Named("Ji",    loglik          ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     status        ),
            Rcpp::Named("backend",    backend       ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", iterations    )
        ))
    );
}

// "VE-step" result: B and Omega fixed, only (M, S2) optimized.
inline Rcpp::List make_vestep_result(
    const arma::mat & M, const arma::mat & S2,
    const arma::vec & loglik,
    int status, const char * backend,
    const std::vector<double> & objective_vec, int iterations
) {
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S2") = S2,
        Rcpp::Named("Ji") = loglik,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     status        ),
            Rcpp::Named("backend",    backend       ),
            Rcpp::Named("objective",  objective_vec ),
            Rcpp::Named("iterations", iterations    )
        ))
    );
}
