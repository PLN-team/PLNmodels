// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

#include "nlopt_wrapper.h"
#include "packer.h"

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
    return -logfact(std::move(y)) + 0.5 * (1. + (1. - double(p)) * std::log(2. * M_PI));
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA_ID>(parameters, init_Theta);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA_ID>(packed, list["Theta"]);
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &Y, &X, &O, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat Theta = packer.unpack<THETA_ID>(parameters);
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) - 0.5 * w_bar * real(log_det(Omega));

        packer.pack<THETA_ID>(grad_storage, (A - Y).t() * (X.each_col() % w));
        packer.pack<M_ID>(grad_storage, diagmat(w) * (M * Omega + A - Y));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat Theta = packer.unpack<THETA_ID>(parameters);
    // Variance parameters
    arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(sum(S2.each_col() % w, 0)));
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::vec>(init_parameters["S"]);         // (n)

    const auto packer = make_packer(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA_ID>(parameters, init_Theta);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA_ID>(packed, list["Theta"]);
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat Theta = packer.unpack<THETA_ID>(parameters);
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::vec S = packer.unpack<S_ID>(parameters);

        arma::vec S2 = S % S;
        const arma::uword p = Y.n_cols;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z.each_col() + 0.5 * S2);
        double sigma2 = arma::as_scalar(accu(M % (M.each_col() % w)) / (w_bar * double(p)) + accu(w % S2) / w_bar);
        double objective = accu(diagmat(w) * (A - Y % Z)) - 0.5 * double(p) * accu(w % log(S2)) +
                           0.5 * w_bar * double(p) * log(sigma2);

        packer.pack<THETA_ID>(grad_storage, (A - Y).t() * (X.each_col() % w));
        packer.pack<M_ID>(grad_storage, diagmat(w) * (M / sigma2 + A - Y));
        packer.pack<S_ID>(grad_storage, w % (S % sum(A, 1) - double(p) * pow(S, -1) - double(p) * S / sigma2));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters); // vec(n) -> mat(n, 1)
    arma::vec S2 = S % S;
    // Regression parameters
    arma::mat Theta = packer.unpack<THETA_ID>(parameters);
    // Variance parameters
    const arma::uword p = Y.n_cols;
    const double n_sigma2 = arma::as_scalar(dot(w, sum(pow(M, 2), 1) + double(p) * S2));
    const double sigma2 = n_sigma2 / (double(p) * w_bar);
    arma::mat Sigma = arma::eye(p, p) * sigma2;
    arma::mat Omega = arma::eye(p, p) * pow(sigma2, -1);
    // Element-wise log-likelihood
    arma::mat Z = O + X * Theta.t() + M;
    arma::mat A = exp(Z.each_col() + 0.5 * S2);
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA_ID>(parameters, init_Theta);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA_ID>(packed, list["Theta"]);
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat Theta = packer.unpack<THETA_ID>(parameters);
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S2);
        arma::rowvec diag_sigma = sum(M % (M.each_col() % w) + (S2.each_col() % w), 0) / w_bar;
        double objective = accu(diagmat(w) * (A - Y % Z - 0.5 * log(S2))) + 0.5 * w_bar * accu(log(diag_sigma));

        packer.pack<THETA_ID>(grad_storage, (A - Y).t() * (X.each_col() % w));
        packer.pack<M_ID>(grad_storage, diagmat(w) * ((M.each_row() / diag_sigma) + A - Y));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S.each_row() % pow(diag_sigma, -1) + S % A - pow(S, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
    arma::mat S2 = S % S;
    // Regression parameters
    arma::mat Theta = packer.unpack<THETA_ID>(parameters);
    // Variance parameters
    arma::rowvec sigma2 = w.t() * (pow(M, 2) + S2) / w_bar;
    arma::vec omega2 = pow(sigma2.t(), -1);
    arma::mat Sigma = diagmat(sigma2);
    arma::mat Omega = diagmat(omega2);
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_B = Rcpp::as<arma::mat>(init_parameters["B"]);         // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,q)

    const auto packer = make_packer(init_Theta, init_B, init_M, init_S);
    enum { THETA_ID, B_ID, M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA_ID>(parameters, init_Theta);
    packer.pack<B_ID>(parameters, init_B);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA_ID>(packed, list["Theta"]);
        packer.pack_double_or_arma<B_ID>(packed, list["B"]);
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat Theta = packer.unpack<THETA_ID>(parameters);
        arma::mat B = packer.unpack<B_ID>(parameters);
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M * B.t();
        arma::mat A = exp(Z + 0.5 * S2 * (B % B).t());
        double objective = accu(diagmat(w) * (A - Y % Z)) + 0.5 * accu(diagmat(w) * (M % M + S2 - log(S2) - 1.));

        packer.pack<THETA_ID>(grad_storage, (A - Y).t() * (X.each_col() % w));
        packer.pack<B_ID>(grad_storage, (diagmat(w) * (A - Y)).t() * M + (A.t() * (S2.each_col() % w)) % B);
        packer.pack<M_ID>(grad_storage, diagmat(w) * ((A - Y) * B + M));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S - 1. / S + A * (B % B) % S));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat Theta = packer.unpack<THETA_ID>(parameters);
    arma::mat B = packer.unpack<B_ID>(parameters);
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
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
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_Theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_M = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_S = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_Theta, init_M, init_S);
    enum { THETA_ID, M_ID, S_ID }; // Names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA_ID>(parameters, init_Theta);
    packer.pack<M_ID>(parameters, init_M);
    packer.pack<S_ID>(parameters, init_S);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA_ID>(packed, list["Theta"]);
        packer.pack_double_or_arma<M_ID>(packed, list["M"]);
        packer.pack_double_or_arma<S_ID>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &O, &X, &Y, &w, &Omega](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat Theta = packer.unpack<THETA_ID>(parameters);
        arma::mat M = packer.unpack<M_ID>(parameters);
        arma::mat S = packer.unpack<S_ID>(parameters);

        arma::mat S2 = S % S;
        arma::mat Z = O + X * Theta.t() + M;
        arma::mat A = exp(Z + 0.5 * S);
        arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
        double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) - trace(Omega * nSigma);

        packer.pack<THETA_ID>(grad_storage, (A - Y).t() * (X.each_col() % w));
        packer.pack<M_ID>(grad_storage, diagmat(w) * (M * Omega + A - Y));
        packer.pack<S_ID>(grad_storage, diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat Theta = packer.unpack<THETA_ID>(parameters);
    arma::mat M = packer.unpack<M_ID>(parameters);
    arma::mat S = packer.unpack<S_ID>(parameters);
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
