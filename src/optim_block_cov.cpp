#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// ---------------------------------------------------------------------------------------
  // Fully parametrized covariance

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_block(
  const Rcpp::List & data  , // List(Y, X, O, w)
  const Rcpp::List & params, // List(B, M, S)
  const Rcpp::List & config  // List of config values
) {
  // Conversion from R, prepare optimization
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const auto init_B   = Rcpp::as<arma::mat>(params["B"]); // (d,p)
  const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  arma::mat Tau       = Rcpp::as<arma::mat>(data["Tau"]); // (q,p)

  const auto metadata = tuple_metadata(init_B, init_M, init_S);
  enum { B_ID, M_ID, S_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<B_ID>(parameters.data()) = init_B;
  metadata.map<M_ID>(parameters.data()) = init_M;
  metadata.map<S_ID>(parameters.data()) = init_S;

  auto optimizer = new_nlopt_optimizer(config, parameters.size());
  if(config.containsElementNamed("xtol_abs")) {
    SEXP value = config["xtol_abs"];
    if(Rcpp::is<double>(value)) {
      set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
    } else {
      auto per_param_list = Rcpp::as<Rcpp::List>(value);
      auto packed = std::vector<double>(metadata.packed_size);
      set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
      set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
      set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
      set_per_value_xtol_abs(optimizer.get(), packed);
    }
  }

  const double w_bar = accu(w);

  // Optimize
  auto objective_and_grad = [&metadata, &Y, &X, &O, &w, &w_bar, &Tau](const double * params, double * grad) -> double {
    const arma::mat B = metadata.map<B_ID>(params);
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    const double w_bar = accu(w);

    arma::mat S2 = S % S;
    arma::mat mu = O + X * B ;
    arma::mat A1 = trunc_exp(M + .5 * S2) ;
    arma::mat A2 = trunc_exp(mu) ;
    arma::mat A = A2 % (A1 * Tau) ;
    arma::mat A_tau = ((A2 * Tau.t()) % A1) ;
    arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    double objective = accu(w.t() * (A - Y % (mu + M * Tau))) - 0.5 * accu(w.t() * log(S2)) - 0.5 * w_bar * real(log_det(Omega));


    metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);
    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A_tau - Y * Tau.t());
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + A_tau % S - pow(S, -1));

    arma::colvec log_alpha = arma::log(mean(Tau,1));
    arma::mat rho = M.t() * Y - A1.t() * A2  ; rho.each_col() += log_alpha ;

    arma::rowvec rho_max = arma::max(rho, 0) ;
    rho_max.print() ;

    rho = arma::trunc_exp(rho.each_row() - rho_max) ;
    Tau = arma::clamp(rho.each_col() / sum(rho, 0), 1e-5, 1-1e-5);

    return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

  // Variational parameters
  arma::mat M = metadata.copy<M_ID>(parameters.data());
  arma::mat S = metadata.copy<S_ID>(parameters.data());
  arma::mat S2 = S % S;
  // Regression parameters
  arma::mat B = metadata.copy<B_ID>(parameters.data());
  // Variance parameters
  arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
  arma::mat Omega = inv_sympd(Sigma);


  // Element-wise log-likehood

  arma::mat mu = O + X * B;
  arma::mat Z = mu + M * Tau;
  arma::mat A1 = trunc_exp(M + .5 * S2) ;
  arma::mat A2 = trunc_exp(mu) ;
  arma::mat A = A2 % (A1 * Tau) ;
  arma::colvec log_alpha = arma::log(mean(Tau,1));
  arma::mat rho = M.t() * Y - A1.t() * A2  ; rho.each_col() += log_alpha ;
  rho = arma::trunc_exp(rho.each_row() - arma::max(rho, 0)) ;
  Tau = arma::clamp(rho.each_col() / sum(rho, 1), 1e-5, 1-1e-5);

  arma::vec loglik = sum(Y % Z - A, 1) + 0.5 * sum(log(S2), 1) - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1) +
  0.5 * real(log_det(Omega)) + ki(Y) + accu(log_alpha.t() * Tau) + accu(Tau % arma::trunc_log(Tau)) ;

  Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
  Ji.attr("weights") = w;
  return Rcpp::List::create(
    Rcpp::Named("B", B),
    Rcpp::Named("M", M),
    Rcpp::Named("S", S),
    Rcpp::Named("Tau", Tau),
    Rcpp::Named("Z", Z),
    Rcpp::Named("A", A),
    Rcpp::Named("Sigma", Sigma),
    Rcpp::Named("Omega", Omega),
    Rcpp::Named("Ji", Ji),
    Rcpp::Named("monitoring", Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("backend", "nlopt"),
        Rcpp::Named("iterations", result.nb_iterations)
    ))
  );
}

// ---------------------------------------------------------------------------------------
  // VE full

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_vestep_block(
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
  const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,p)
  const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,p)

  const auto metadata = tuple_metadata(init_M, init_S);
  enum { M_ID, S_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<M_ID>(parameters.data()) = init_M;
  metadata.map<S_ID>(parameters.data()) = init_S;

  auto optimizer = new_nlopt_optimizer(config, parameters.size());
  if(config.containsElementNamed("xtol_abs")) {
    SEXP value = config["xtol_abs"];
    if(Rcpp::is<double>(value)) {
      set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
    } else {
      auto per_param_list = Rcpp::as<Rcpp::List>(value);
      auto packed = std::vector<double>(metadata.packed_size);
      set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
      set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
      set_per_value_xtol_abs(optimizer.get(), packed);
    }
  }

  // Optimize
  auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &B, &Omega](const double * params, double * grad) -> double {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);

    arma::mat S2 = S % S;
    arma::mat Z = O + X * B + M;
    arma::mat A = exp(Z + 0.5 * S2);
    arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2) ;
    double objective = accu(w.t() * (A - Y % Z - 0.5 * log(S2))) + 0.5 * trace(Omega * nSigma) ;

    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A - Y);
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A - pow(S, -1));

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
    Rcpp::Named("iterations", result.nb_iterations)
  ))
);
}