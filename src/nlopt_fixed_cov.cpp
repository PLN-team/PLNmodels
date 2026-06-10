#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "nlopt_impl.h"

// ---------------------------------------------------------------------------------------
// Fixed covariance PLN — nlopt/CCSAQ optimizer: B profiled via closed form, reduced parameter vector

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_fixed(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    const auto init_B = Rcpp::as<arma::mat>(params["B"]);
    const auto init_M = Rcpp::as<arma::mat>(params["M"]);
    const auto  Omega = Rcpp::as<arma::mat>(params["Omega"]);
    const auto init_S = Rcpp::as<arma::mat>(params["S"]);

    const auto metadata = tuple_metadata(init_M, init_S);
    enum { M_ID, S_ID };

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = arma::log(init_S % init_S);

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    std::vector<double> objective_vec;

    const arma::mat Xw         = X.each_col() % w;
    const arma::mat P_X        = (X.n_cols > 0) ? arma::solve(X.t() * Xw, Xw.t()) : arma::mat(0, Y.n_rows);
    const arma::vec Omega_diag = diagvec(Omega);

    auto objective_and_grad = [&](const double * par, double * grad) -> double {
        const arma::mat M_full = metadata.map<M_ID>(par);
        const arma::mat logS2  = metadata.map<S_ID>(par);
        const arma::mat B      = P_X * M_full;
        const arma::mat M_res  = M_full - X * B;
        arma::mat gM, gS;
        const double obj = full_cov_obj_grad_impl(M_res, O + M_full, logS2, Omega, Omega_diag, Y, w, gM, gS);
        metadata.map<M_ID>(grad) = gM;
        metadata.map<S_ID>(grad) = gS;
        objective_vec.push_back(obj);
        return obj;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M      = metadata.copy<M_ID>(parameters.data());  // M_full
    arma::mat logS2  = metadata.copy<S_ID>(parameters.data());
    arma::mat S2     = arma::exp(logS2);
    arma::mat S      = arma::exp(0.5 * logS2);
    arma::mat B      = P_X * M;
    arma::mat M_res  = M - X * B;
    arma::mat Sigma  = (M_res.t() * (M_res.each_col() % w) + diagmat(w.t() * S2)) / accu(w);
    arma::mat Z      = O + M;
    arma::mat A      = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * ((M_res * Omega) % M_res - logS2 + S2 * diagmat(Omega)), 1)
                     + 0.5 * real(log_det(Omega)) + ki(Y);

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
