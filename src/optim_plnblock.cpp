#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// [[Rcpp::export]]
arma::vec plnblock_vloglik(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat &T  = Rcpp::as<arma::mat>(params["T"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword n = Y.n_rows;
  const arma::uword p = Y.n_cols;
  const arma::uword q = M.n_cols;
  const arma::mat S2 = S % S ;
  const arma::mat mu = O + X * B;
  const arma::mat Z = mu + M * T;
  const arma::mat A = trunc_exp(mu) % (trunc_exp(M + .5 * S2) * T) ;

  // Element-wise log-likehood
  return(
    0.5 * real(log_det(Omega)) + 0.5 * double(q) - logfact(Y)
  + sum(Y % Z - A, 1) + 0.5 * sum(log(S2), 1)
  - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)
   + accu(log_pi.t() * T) - accu(T % arma::trunc_log(T))
  ) ;
}

// [[Rcpp::export]]
arma::vec plnblock_loglik(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat &T  = Rcpp::as<arma::mat>(params["T"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword n = Y.n_rows;
  const arma::uword p = Y.n_cols;
  const arma::uword q = M.n_cols;
  const arma::mat S2 = S % S ;
  const arma::mat mu = O + X * B;
  const arma::mat Z = mu + M * T;
  const arma::mat A = trunc_exp(mu) % (trunc_exp(M + .5 * S2) * T) ;

  // Element-wise log-likehood
  return(
   (w.t() * (
    0.5 * real(log_det(Omega)) + 0.5 * double(q) - logfact(Y)
    + sum(Y % Z - A, 1) + 0.5 * sum(log(S2), 1)
    - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)
   )) + accu(log_pi.t() * T) - accu(T % arma::trunc_log(T))
  ) ;
}

// [[Rcpp::export]]
Rcpp::List  optim_plnblock_Omega(
  const arma::mat & M, // (n,q)
  const arma::mat & S, // (n,q)
  const arma::vec & w  // (n)
) {
  double w_bar = accu(w);
  arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * (S % S))) ;

  return Rcpp::List::create(
    Rcpp::Named("Sigma") = (1./w_bar) * nSigma,
    Rcpp::Named("Omega") = w_bar * inv_sympd(nSigma));
}

// [[Rcpp::export]]
arma::mat optim_plnblock_Tau(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));
  const arma::mat mu = O + X * B ;
  const arma::mat A1 = trunc_exp(M + .5 * S % S ) ;
  const arma::mat A2 = trunc_exp(mu) ;

  arma::mat Tau = M.t() * Y - A1.t() * A2  ;
  Tau.each_col() += log_pi ;
  Tau.each_col( [](arma::vec& x){
    x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
  }) ;
  return Tau ;
}

// [[Rcpp::export]]
Rcpp::List optim_plnblock_B(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega)
    const Rcpp::List & configuration // List of config values
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const auto init_B   = Rcpp::as<arma::mat>(params["B"]); // (n,p)
  const auto metadata = tuple_metadata(init_B);
  enum { B_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<B_ID>(parameters.data()) = init_B;
  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
  set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));

  const arma::mat A1_T = trunc_exp(M + .5 * S % S) * T;

  // Optimize
  auto objective_and_grad =
    [&metadata, &Y, &X, &O, &A1_T, &T, &w](const double * params, double * grad) -> double {
      const arma::mat B = metadata.map<B_ID>(params);

      arma::mat mu = O + X * B ;
      arma::mat A = trunc_exp(mu) % A1_T ;

      double objective = accu(w.t() * (A - Y % mu));

      metadata.map<B_ID>(grad) = (X.each_col() % w).t() * (A - Y);

      return objective;

    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat B = metadata.copy<B_ID>(parameters.data());
    return Rcpp::List::create(
      Rcpp::Named("status") = static_cast<int>(result.status),
      Rcpp::Named("iterations") = result.nb_iterations,
      Rcpp::Named("objective") = result.objective,
      Rcpp::Named("B") = B);
}

// [[Rcpp::export]]
Rcpp::List optim_plnblock_VE(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega)
    const Rcpp::List & configuration // List of config values
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const auto metadata = tuple_metadata(init_M, init_S);
  enum { M_ID, S_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<M_ID>(parameters.data()) = init_M;
  metadata.map<S_ID>(parameters.data()) = init_S;
  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
  set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));

  const arma::mat mu = O + X * B ;
  const arma::mat A2 = trunc_exp(mu) ;

  // Optimize
  auto objective_and_grad = [&metadata, &Y, &X, &mu, &A2,  &T, &Omega, &w](const double * params, double * grad) -> double {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);

    arma::mat S2 = S % S ;
    arma::mat A1  = trunc_exp(M + .5 * S2) ;
    arma::mat A   = A2 % (A1 * T) ;
    arma::mat A_T = A1 % (A2 * T.t()) ;
    arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
    double objective =
      accu(w.t() * (A - Y % (mu + M * T))) - 0.5 * accu(w.t() * log(S2))
      + .5 * trace(Omega * nSigma) ;

    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A_T - Y * T.t());
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A_T - pow(S, -1));

    return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    return Rcpp::List::create(
      Rcpp::Named("status") = static_cast<int>(result.status),
      Rcpp::Named("iterations") = result.nb_iterations,
      Rcpp::Named("objective") = result.objective,
      Rcpp::Named("M") = M,
      Rcpp::Named("S") = S);
}
