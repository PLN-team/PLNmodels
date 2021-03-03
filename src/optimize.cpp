// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

#include "nlopt_wrapper.h"
#include "packing.h"

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
//    return -logfact(std::move(y)) + 0.5 * (1. + (1. - double(p)) * std::log(2. * M_PI));
    return -logfact(std::move(y)) + 0.5 * double(p) ;
}

// ---------------------------------------------------------------------------------------
// Fully parametrized covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_full(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

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
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad = [&metadata, &Y, &X, &O, &w, &w_bar](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) - 0.5 * w_bar * real(log_det(Omega));

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));

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
    arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    arma::mat Omega = inv_sympd(Sigma);
    // Element-wise log-likehood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::vec loglik = sum(Y % Z - A + 0.5 * log(S2) - 0.5 * ((M * Omega) % M + S2 * diagmat(Omega)), 1) +
                       0.5 * real(log_det(Omega)) + ki(Y);

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
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Spherical covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_spherical(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::vec>(init_parameters["S"]);         // (n)

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

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
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &w_bar](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2 * arma::ones(p).t());
        double sigma2 = arma::as_scalar(dot(w, sum(pow(M, 2), 1) + double(p) * S2) / (double(p) * w_bar) );
        double objective = accu(w.t() * (A - Y % Z)) - 0.5 * double(p) * dot(w, log(S2/sigma2)) ;

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * (M / sigma2 + A - Y);
        metadata.map<S_ID>(grad) = w % (double(p) * S / sigma2 + S % sum(A, 1) - double(p) * pow(S, -1));

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Variational parameters
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data()); // vec(n) -> mat(n, 1)
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    // Variance parameters
    const arma::uword p = Y.n_cols;
    const double sigma2 = arma::as_scalar(dot(w, sum(pow(M, 2), 1) + double(p) * S2)) / (double(p) * w_bar);
    arma::sp_mat Sigma(p,p); Sigma.diag() = arma::ones<arma::vec>(p) * sigma2;
    arma::sp_mat Omega(p,p); Omega.diag() = arma::ones<arma::vec>(p) * pow(sigma2, -1);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2 * arma::ones(p).t());
    arma::mat loglik = sum(Y % Z - A - 0.5 * pow(M, 2) / sigma2, 1) - 0.5 * double(p) * S2 / sigma2 +
                       0.5 * double(p) * log(S2 / sigma2) + ki(Y);

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
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Diagonal covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_diagonal(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

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
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &w_bar](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::rowvec diag_sigma = w.t() * (M % M + S2) / w_bar;
        double objective = accu(diagmat(w) * (A - Y % Z - 0.5 * log(S2))) + 0.5 * w_bar * accu(log(diag_sigma));

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * ((M.each_row() / diag_sigma) + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % pow(diag_sigma, -1) + S % A - pow(S, -1));
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
    arma::rowvec sigma2 = w.t() * (M % M + S2) / w_bar;
    arma::vec omega2 = pow(sigma2.t(), -1);
    arma::sp_mat Sigma(Y.n_cols, Y.n_cols);
    Sigma.diag() = sigma2;
    arma::sp_mat Omega(Y.n_cols, Y.n_cols);
    Omega.diag() = omega2;
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik =
        sum(Y % Z - A + 0.5 * log(S2), 1) - 0.5 * (pow(M, 2) + S2) * omega2 + 0.5 * sum(log(omega2)) + ki(Y);

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
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List cpp_optimize_rank(
    const Rcpp::List & init_parameters, // List(Theta, B, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_B = Rcpp::as<arma::mat>(init_parameters["B"]);         // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,q)

    const auto metadata = tuple_metadata(init_Theta, init_B, init_M, init_S);
    enum { THETA_ID, B_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
    if(configuration.containsElementNamed("xtol_abs")) {
        SEXP value = configuration["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<THETA_ID>(packed.data()), per_param_list["Theta"]);
            set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M * B.t();
        arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
        double objective = accu(diagmat(w) * (A - Y % Z)) + 0.5 * accu(diagmat(w) * (M % M + S2 - log(S2) - 1.));

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<B_ID>(grad) = (diagmat(w) * (A - Y)).t() * M + (A.t() * (S2.each_col() % w)) % B;
        metadata.map<M_ID>(grad) = diagmat(w) * ((A - Y) * B + M);
        metadata.map<S_ID>(grad) = diagmat(w) * (S - 1. / S + A * (B % B) % S);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat Sigma = B * (M.t() * (M.each_col() % w) + diagmat(sum(S2.each_col() % w, 0))) * B.t() / accu(w);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M * B.t();
    arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
    arma::mat loglik = arma::sum(Y % Z - A, 1) - 0.5 * sum(M % M + S2 - log(S2) - 1., 1) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", Theta),
        Rcpp::Named("B", B),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Sparse inverse covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_sparse(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & Y,                // responses (n,p)
    const arma::mat & X,                // covariates (n,d)
    const arma::mat & O,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const arma::mat & Omega,            // covinv (p,p)
    const Rcpp::List & configuration    // List of config values
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto metadata = tuple_metadata(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<THETA_ID>(parameters.data()) = init_Theta;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

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
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &Omega](const double * params, double * grad) -> double {
        const arma::mat Theta = metadata.map<THETA_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) - trace(Omega * nSigma);

        metadata.map<THETA_ID>(grad) = (A - Y).t() * (X.each_col() % w);
        metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
        metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat Theta = metadata.copy<THETA_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat Sigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) / accu(w);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat loglik = sum(Y % Z - A - 0.5 * ((M * Omega) % M - log(S2) + S2 * diagmat(Omega)), 1) +
                       0.5 * real(log_det(Omega)) + ki(Y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", Theta),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("loglik", loglik));
}
