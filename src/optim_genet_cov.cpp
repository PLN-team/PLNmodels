#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
// Covariance for heritability

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_genetic_modeling(
    const Rcpp::List & init_parameters, // List(Theta, M, S, rho)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const arma::mat & C,                // correlation matrix (p x p)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M     = Rcpp::as<arma::mat>(init_parameters["M"]);     // (n,p)
    const auto init_S     = Rcpp::as<arma::mat>(init_parameters["S"]);     // (n,p)
    const auto init_rho   = Rcpp::as<double>(init_parameters["rho"]);      // double

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S, init_rho);
    enum { THETA_ID, M_ID, S_ID, RHO_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;
    metadata.map<RHO_ID>(parameters.data()) = init_rho;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<THETA_ID>(packed.data()), per_param_list["Theta"]);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_from_r_sexp(metadata.map<RHO_ID>(packed.data()), per_param_list["rho"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    // Some fixed quantities along optimization
    const double w_bar = accu(w);
    // Diagonalization of matrix C
    arma::vec Lambda;
    arma::mat V;
    arma::eig_sym(Lambda, V, C);

    // Optimize
    auto objective_and_grad = [&metadata, &Y, &X, &O, &V, &Lambda, &w, &w_bar](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        const double rho  = metadata.map<RHO_ID>(params);

        const arma::uword p = Y.n_cols;
        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::vec u = rho * Lambda + (1-rho) ;
        arma::mat R = V.t() * (M.t() * M + diagmat(w.t() * S2)) * V ;
        double sigma2 = accu( diagvec(R) / u ) / (double(p) * w_bar);
        arma::mat Omega = V * diagmat(pow(sigma2 * u, -1)) * V.t() ;

        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) +
              0.5 * trace(Omega * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2))) +
              0.5 * w_bar * accu(log(u * sigma2));

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));
        metadata.map<RHO_ID>(grad) = accu(0.5 * w_bar * (Lambda - 1) / u - (0.5/sigma2) * diagvec(R) % (Lambda - 1) / pow(u, 2) );

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    // Variance parameters
    const arma::uword p = Y.n_cols;
    double rho = metadata.copy<RHO_ID>(parameters.data()) ;

    arma::vec u = rho * Lambda + (1-rho) ;
    arma::mat R = V.t() * (M.t() * M + diagmat(w.t() * S2)) * V ;
    double sigma2 = accu( diagvec(R) / u ) / (double(p) * w_bar);
    arma::mat Omega = V * diagmat(pow(sigma2 * u, -1)) * V.t() ;
    arma::mat Sigma = V * diagmat(sigma2 * u) * V.t() ;
    // Element-wise log-likehood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * log(S2) - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) -
                       0.5 * accu(log(sigma2*u)) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", Theta),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("sigma2", sigma2),
        Rcpp::Named("rho", rho),
        Rcpp::Named("loglik", loglik));
}

