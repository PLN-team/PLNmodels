#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Rank-constrained covariance
//
// Variational parameter: ψ = log(S²) (unconstrained) instead of S (bounded > 0).
// Interface: takes S2 (variance), converts to ψ = log(S²) before optimizing, returns S2.
// Gradient w.r.t. ψ: ½ · w ⊙ (S² ⊙ (1 + A·C²) − 1)   [no 1/S singularity]

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_rank(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, C, M, S2)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const PlnData d(data);
    const auto init_B   = Rcpp::as<arma::mat>(params["B"]);  // (d,p)
    const auto init_C   = Rcpp::as<arma::mat>(params["C"]);  // (p,q)
    const auto init_M   = Rcpp::as<arma::mat>(params["M"]);  // (n,q)
    const arma::mat init_S2  = Rcpp::as<arma::mat>(params["S2"]);       // (n,q) variance
    const arma::mat init_psi = arma::log(init_S2);                      // ψ = log(S²)

    const auto metadata = tuple_metadata(init_B, init_C, init_M, init_psi);
    enum { B_ID, C_ID, M_ID, PSI_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>  (parameters.data()) = init_B;
    metadata.map<C_ID>  (parameters.data()) = init_C;
    metadata.map<M_ID>  (parameters.data()) = init_M;
    metadata.map<PSI_ID>(parameters.data()) = init_psi;

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));

    const arma::mat Xw = d.X.each_col() % d.w;   // precomputed once

    auto objective_and_grad = [&metadata, &d, &Xw, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat B   = metadata.map<B_ID>  (params);
        const arma::mat C   = metadata.map<C_ID>  (params);
        const arma::mat M   = metadata.map<M_ID>  (params);
        const arma::mat psi = metadata.map<PSI_ID>(params);

        const arma::mat C2 = C % C;
        const arma::mat S2 = arma::exp(psi);                             // S² = exp(ψ)
        arma::mat Z = d.O + d.X * B + M * C.t();
        arma::mat A = arma::exp(Z + 0.5 * S2 * C2.t());
        const arma::mat AmY = A - d.Y;  // gradient uses (A-Y); objective uses (A - Y%Z)
        double objective = accu(d.w.t() * arma::sum(A - d.Y % Z, 1))
                         + 0.5 * accu(d.w.t() * arma::sum(M % M + S2 - psi - 1., 1));

        arma::mat gC  = (AmY.each_col() % d.w).t() * M + (A.t() * (S2.each_col() % d.w)) % C;
        arma::mat gM  = AmY * C + M;  gM.each_col()  %= d.w;
        arma::mat gPS = 0.5 * (S2 % (1. + A * C2) - 1.); gPS.each_col() %= d.w;
        metadata.map<B_ID>  (grad) = Xw.t() * AmY;
        metadata.map<C_ID>  (grad) = gC;
        metadata.map<M_ID>  (grad) = gM;
        metadata.map<PSI_ID>(grad) = gPS;

        objective_vec.push_back(objective);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B   = metadata.copy<B_ID>  (parameters.data());
    arma::mat C   = metadata.copy<C_ID>  (parameters.data());
    arma::mat M   = metadata.copy<M_ID>  (parameters.data());
    arma::mat psi = metadata.copy<PSI_ID>(parameters.data());
    arma::mat S2  = arma::exp(psi);
    arma::mat Sigma = C * (M.t() * (M.each_col() % d.w) + arma::diagmat(arma::sum(S2.each_col() % d.w, 0))) * C.t() / accu(d.w);
    arma::mat Omega = C * inv_sympd((M.t() * (M.each_col() % d.w) + arma::diagmat(arma::sum(S2.each_col() % d.w, 0))) / accu(d.w)) * C.t();
    arma::mat Z   = d.O + d.X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * (C % C).t());
    arma::mat loglik = arma::sum(d.Y % Z - A, 1) - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(d.Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = d.w;
    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S2", S2),
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

// ---------------------------------------------------------------------------------------
// VE rank
// Rank-constrained covariance (for prediction in the PCA space)

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_rank(
        const Rcpp::List & data  , // List(Y, X, O, w)
        const Rcpp::List & params, // List(M, S2)
        const arma::mat & B,       // (d,p)
        const arma::mat & C,       // (p,q)
        const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const PlnData d(data);
    const auto init_M        = Rcpp::as<arma::mat>(params["M"]);  // (n,q)
    const arma::mat init_S2  = Rcpp::as<arma::mat>(params["S2"]); // (n,q) variance
    const arma::mat init_psi = arma::log(init_S2);                // ψ = log(S²)

    const auto metadata = tuple_metadata(init_M, init_psi);
    enum { M_ID, PSI_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>  (parameters.data()) = init_M;
    metadata.map<PSI_ID>(parameters.data()) = init_psi;

    // Optimize
    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;
    objective_vec.reserve(nlopt_get_maxeval(optimizer.get()));
    const arma::mat C2 = C % C;
    const arma::mat XB = d.X * B;

    auto objective_and_grad = [&metadata, &d, &XB, &C, &C2, &objective_vec](const double * params, double * grad) -> double {
        const arma::mat M   = metadata.map<M_ID>  (params);
        const arma::mat psi = metadata.map<PSI_ID>(params);

        arma::mat S2 = arma::exp(psi);
        arma::mat Z  = d.O + XB + M * C.t();
        arma::mat A  = arma::exp(Z + 0.5 * S2 * C2.t());
        const arma::mat AmY = A - d.Y;  // gradient uses (A-Y); objective uses (A - Y%Z)
        double objective = accu(d.w.t() * arma::sum(A - d.Y % Z, 1))
                         + 0.5 * accu(d.w.t() * arma::sum(M % M + S2 - psi - 1., 1));

        arma::mat gM  = AmY * C + M;  gM.each_col()  %= d.w;
        arma::mat gPS = 0.5 * (S2 % (1. + A * C2) - 1.); gPS.each_col() %= d.w;
        metadata.map<M_ID>  (grad) = gM;
        metadata.map<PSI_ID>(grad) = gPS;

        objective_vec.push_back(objective);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat M   = metadata.copy<M_ID>  (parameters.data());
    arma::mat psi = metadata.copy<PSI_ID>(parameters.data());
    arma::mat S2  = arma::exp(psi);
    arma::mat Z   = d.O + d.X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    arma::mat loglik = arma::sum(d.Y % Z - A, 1) - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(d.Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = d.w;
    return make_vestep_result(M, S2, Ji, static_cast<int>(result.status), "nlopt", objective_vec, result.nb_iterations);
}
