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
  const double w_bar  = accu(w);
  const auto init_B   = Rcpp::as<arma::mat>(params["B"]); // (d,p)
  const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  arma::mat Tau = Rcpp::as<arma::mat>(params["Tau"]);     // (q,p)

  const auto metadata = tuple_metadata(init_B, init_M, init_S);
  enum { B_ID, M_ID, S_ID}; // Names for metadata indexes

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

  // Optimization over B, M and S: computation of objective and gradients
  auto objective_and_grad = [&metadata, &Y, &X, &O, &w, &w_bar, &Tau](const double * params, double * grad) -> double {
    const arma::mat B = metadata.map<B_ID>(params);
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);

    arma::mat S2 = S % S;
    arma::mat mu = O + X * B ;
    arma::mat A1 = trunc_exp(M + .5 * S2) ;
    arma::mat A2 = trunc_exp(mu) ;

    arma::mat A = A2 % (A1 * Tau) ;
    arma::mat A_tau = A1 % (A2 * Tau.t()) ;
    arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    double objective = accu(w.t() * (A - Y % (mu + M * Tau))) - 0.5 * accu(w.t() * log(S2)) - 0.5 * w_bar * real(log_det(Omega));

    metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);
    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A_tau - Y * Tau.t());
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A_tau - pow(S, -1));

    return objective;
  };

  // Posterior probabilities (Tau): compute updates
  // arma::vec param
  auto update_posterior_prob = [&metadata, &Y, &X, &O, &Tau](const double *params) {
    const arma::mat B = metadata.map<B_ID>(params);
    const arma::vec log_alpha = arma::log(mean(Tau,1));
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    arma::mat S2 = S % S;
    arma::mat mu = O + X * B ;
    arma::mat A1 = trunc_exp(M + .5 * S2) ;
    arma::mat A2 = trunc_exp(mu) ;

    Tau = M.t() * Y - A1.t() * A2  ;
    Tau.each_col() += log_alpha ;
    Tau.each_col( [](arma::vec& x){
      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
    }) ;
  };

  // =========================================================================
  //
  // MAIN OPTIMIZATION

  // Initialization and first step
  int iter = 0;
  Rcpp::List posteriorProb ; posteriorProb.push_back(arma::mat(Tau)) ;
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
  Rcpp::NumericVector objective ; objective.push_back(result.objective) ;
  // std::cout << "objective" << result.objective << "  inner iteration" << result.nb_iterations << std::endl;

  update_posterior_prob(parameters.data()) ;

  // Alternate optimization
  double threshold = Rcpp::as<double>(config["ftol_out"]);
  double maxiter = Rcpp::as<int>(config["maxit_out"]);
  do { // GO FOR IT!!
    iter++;

    // Optimize Tau
    update_posterior_prob(parameters.data()) ;

    // Optimize B, M and S
    result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
    objective.push_back(result.objective)  ;

    // std::cout << "objective" << result.objective << "  inner iteration" << result.nb_iterations << std::endl;

  } while( std::abs(objective[iter] - objective[iter-1])/std::abs(objective[iter-1]) > threshold & iter+1 < maxiter) ;

  // =========================================================================
  //
  // MAIN OPTIMIZATION DONE
  //
  // Preparing output

  // Variational parameters
  arma::mat M = metadata.copy<M_ID>(parameters.data());
  arma::mat S = metadata.copy<S_ID>(parameters.data());
  arma::mat S2 = S % S;
  // Regression parameters
  arma::mat B = metadata.copy<B_ID>(parameters.data());
  // Variance parameters
  arma::mat Sigma = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
  arma::mat Omega = inv_sympd(Sigma);
  // Group proportion
  arma::vec log_alpha = arma::log(mean(Tau, 1));
  // Latent variable and expected value
  arma::mat mu = O + X * B;
  arma::mat Z = mu + M * Tau;
  arma::mat A = trunc_exp(mu) % (trunc_exp(M + .5 * S2) * Tau) ;

  // Element-wise log-likehood
  arma::vec loglik = sum(Y % Z - A, 1) + 0.5 * sum(log(S2), 1) - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1) +
  0.5 * real(log_det(Omega)) - logfact(Y) + .5 * M.n_cols + accu(log_alpha.t() * Tau) - accu(Tau % arma::trunc_log(Tau)) ;

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
        // Rcpp::Named("posteriorProb", wrap(posteriorProb)),
        Rcpp::Named("objective", wrap(objective)),
        Rcpp::Named("outer_iterations", objective.length()),
        Rcpp::Named("iterations", result.nb_iterations)
    ))
  );
}

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_block_sparse(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, M, S)
    const Rcpp::List & config  // List of config values
) {
  // Conversion from R, prepare optimization
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const auto init_B   = Rcpp::as<arma::mat>(params["B"]);   // (d,p)
  const auto init_M   = Rcpp::as<arma::mat>(params["M"]);   // (n,q)
  const auto init_S   = Rcpp::as<arma::mat>(params["S"]);   // (n,q)
  arma::mat Tau       = Rcpp::as<arma::mat>(params["Tau"]); // (q,p)
  const Rcpp::NumericMatrix rho = params["rho"] ; // double
  const double w_bar = accu(w);
  arma::mat Omega ; // (q,q)
  arma::mat Sigma ; // (q,q)

  const auto metadata = tuple_metadata(init_B, init_M, init_S);
  enum { B_ID, M_ID, S_ID}; // Names for metadata indexes

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

  // Optimization over B, M and S: computation of objective and gradients (with Omega fixed)
  auto objective_and_grad = [&metadata, &Y, &X, &O, &w, &w_bar, &Omega, &Tau](const double * params, double * grad) -> double {
    const arma::mat B = metadata.map<B_ID>(params);
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    const double w_bar = accu(w);

    arma::mat S2 = S % S;
    arma::mat mu = O + X * B ;
    arma::mat A1 = trunc_exp(M + .5 * S2) ;
    arma::mat A2 = trunc_exp(mu) ;

    arma::mat A = A2 % (A1 * Tau) ;
    arma::mat A_tau = A1 % (A2 * Tau.t()) ;
    arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
    double objective = accu(w.t() * (A - Y % (mu + M * Tau))) - 0.5 * accu(w.t() * log(S2)) + 0.5 * trace(Omega * nSigma);

    metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);
    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A_tau - Y * Tau.t());
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A_tau - pow(S, -1));

    return objective;
  };

  // Posterior probabilities (Tau): compute updates
  auto update_posterior_prob = [&metadata, &Y, &X, &O, &Tau](const double * params) -> void {
    const arma::mat B = metadata.map<B_ID>(params);
    const arma::vec log_alpha = arma::log(mean(Tau,1));
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    arma::mat S2 = S % S;
    arma::mat mu = O + X * B ;
    arma::mat A1 = trunc_exp(M + .5 * S2) ;
    arma::mat A2 = trunc_exp(mu) ;

    Tau = M.t() * Y - A1.t() * A2  ;
    Tau.each_col() += log_alpha ;
    Tau.each_col( [](arma::vec& x){
      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
    }) ;
  };

  // Posterior probabilities (Tau): compute updates
  auto graphical_lasso = [&metadata, &rho, &w, &w_bar, &Omega, &Sigma](const double * params) -> void {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    arma::mat S2 = S % S;
    arma::mat Sigma_hat = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
    // Graphical-Lasso: interface to glassoFast function
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("glassoFast");
    Rcpp::Function glassoFast = pkg["glassoFast"];
    Rcpp::List out_glasso = glassoFast(Sigma_hat, rho) ;
    Omega = Rcpp::as<arma::mat>(out_glasso["wi"]) ;
    Sigma = Rcpp::as<arma::mat>(out_glasso["w"]) ;
  };

  // =========================================================================
  //
  // MAIN OPTIMIZATION

  // Initialization and first step
  int iter = 0;
  Rcpp::List posteriorProb ; posteriorProb.push_back(arma::mat(Tau)) ;

  // Omega (first pass of graphical-lasso)
  graphical_lasso(parameters.data()) ;
  // B, M, S: gradient ascent along the ELBO
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
  Rcpp::NumericVector objective ; objective.push_back(result.objective) ;

  // Alternate optimization
  double threshold = Rcpp::as<double>(config["ftol_out"]);
  double maxiter = Rcpp::as<int>(config["maxit_out"]);
  do { // GO for it
    iter++;

    // Optimize Tau
    auto update_posterior_prob(parameters) ;
    posteriorProb.push_back(arma::mat(Tau)) ;

    // Omega (first pass of graphical-lasso)
    auto graphical_lasso(parameters) ;

    // Optimize B, M and S
    result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
    objective.push_back(result.objective)  ;

  } while( std::abs(objective[iter] - objective[iter-1])/std::abs(objective[iter-1]) > threshold & iter < maxiter) ;


  // =========================================================================
  //
  // MAIN OPTIMIZATION DONE
  //
  // Preparing output

  // Variational parameters
  arma::mat M = metadata.copy<M_ID>(parameters.data());
  arma::mat S = metadata.copy<S_ID>(parameters.data());
  arma::mat S2 = S % S;
  // Regression parameters
  arma::mat B = metadata.copy<B_ID>(parameters.data());
  // Variance parameters
  arma::mat Sigma_hat = (1. / w_bar) * (M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
  // Group proportion
  arma::vec log_alpha = arma::log(mean(Tau, 1));
  // Latent variable and expected value
  arma::mat mu = O + X * B;
  arma::mat Z = mu + M * Tau;
  arma::mat A = trunc_exp(mu) % (trunc_exp(M + .5 * S2) * Tau) ;

  // Element-wise log-likehood
  arma::vec loglik = sum(Y % Z - A, 1) + 0.5 * sum(log(S2), 1) - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1) +
    0.5 * real(log_det(Omega)) - logfact(Y) + .5 * M.n_cols + accu(log_alpha.t() * Tau) - accu(Tau % arma::trunc_log(Tau)) ;

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
        Rcpp::Named("posteriorProb", wrap(posteriorProb)),
        Rcpp::Named("objective", wrap(objective)),
        Rcpp::Named("outer_iterations", objective.length()),
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
