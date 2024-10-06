#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

// [[Rcpp::export]]
arma::vec plnblockbis_vloglik(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["Tau"]); // (q,p)
  const arma::rowvec & dm = Rcpp::as<arma::rowvec>(params["dm"]); // (p,p)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)

  const double p = double(Y.n_cols);
  const double q = double(M.n_cols);
  const arma::vec log_pi = arma::trunc_log(mean(T,1));
  const arma::mat S2 = S % S ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat Z = O + Mu + M * T;
  const arma::mat A = trunc_exp(O + Mu + .5 * Delta2) % (trunc_exp(M + .5 * S2) * T) ;

  // Element-wise log-likelihood
  return(
         sum(Y % Z - A, 1) - logfact(Y)                      // Poisson term
      + 0.5 * (real(log_det(Omega)) + sum(log(S2), 1) + q)   // Gaussian-block term
      - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1) // ...
      + 0.5 * (sum(log(dm)) + sum(log(Delta2), 1) + p)       // Gaussian-species term
      - 0.5 * (pow(Mu - X * B, 2) + Delta2) * dm.t()         // ...
      + accu(log_pi.t() * T) - accu(T % arma::trunc_log(T))  // Multinomial-block term
  ) ;
}

// [[Rcpp::export]]
double plnblockbis_loglik(
    const Rcpp::List & data,  // List(Y, X, O, w)
    const Rcpp::List & params // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & T = Rcpp::as<arma::mat>(params["Tau"]); // (n,p)
  const arma::rowvec & dm = Rcpp::as<arma::rowvec>(params["dm"]); // (p,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)

  const double n = double(Y.n_rows);
  const double p = double(Y.n_cols);
  const double q = double(M.n_cols);
  const arma::vec log_pi = arma::trunc_log(mean(T,1));
  const arma::mat S2 = S % S ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat Z = O + Mu + M * T;
  const arma::mat A = trunc_exp(O + Mu + .5 * Delta2) % (trunc_exp(M + .5 * S2) * T) ;
  const arma::mat Sigma = (M.t() * M + diagmat(sum(S2,0))) / n;

  // At optimum, accu(d % dm) = p so some terms vanish
  return (
      accu(Y % (O + Mu + M * T) - A) - accu(logfact(Y))     // Poisson term
    - .5 * n * (accu(Omega % Sigma) - log_det_sympd(Omega)) // Gaussian-block term
    + .5 * n * (accu(log(dm)))                              // Gaussian-species term
    + .5 * (accu(log(S2)) + accu(log(Delta2)) + n * q)      // Gaussian entropies + cst.
    + (accu(log_pi.t() * T) - accu(T % trunc_log(T)))       // multinomial term
  );
}

// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_Omega(
  const arma::mat & M, // (n,q)
  const arma::mat & S  // (n,q)
) {
  arma::mat nSigma = M.t() * M + diagmat(sum(S % S, 0)) ;
  const double n = double(M.n_rows);
  return Rcpp::List::create(
    Rcpp::Named("Sigma") = (1./n) * nSigma,
    Rcpp::Named("Omega") = n * inv_sympd(nSigma));
}

// [[Rcpp::export]]
arma::mat  optim_plnblockbis_B(
    const arma::mat & XtXm,
    const arma::mat & X,
    const arma::mat & Mu
) {
  return(XtXm * X.t() * Mu);
}

// [[Rcpp::export]]
arma::rowvec  optim_plnblockbis_dm(
    const arma::mat & X,
    const arma::mat & B,
    const arma::mat & Mu,
    const arma::mat & Delta
) {
  return(pow(mean(pow(Mu - X * B, 2) + Delta % Delta, 0), -1)) ;
}

// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_Tau(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"])  ; // responses (n,p)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"])  ; // offsets (n,p)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & Tau0 = Rcpp::as<arma::mat>(params["Tau"]); // (n,p)
  const arma::vec log_pi = arma::trunc_log(mean(Tau0,1));
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat A1 = trunc_exp(O + Mu + .5* Delta % Delta) ;
  const arma::mat A2 = trunc_exp(M + .5 * S % S) ;
  const double n = double(Y.n_rows) ;
  arma::mat Tau = (M.t() * Y - A2.t() * A1) / n ;
  Tau.each_col() += log_pi - 1 ;
  Tau.each_col( [](arma::vec& x){
    x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
  }) ;
  return Rcpp::List::create(
    Rcpp::Named("Tau") = Tau,
    Rcpp::Named("A") = A1 % (A2 * Tau)
  ) ;
}

////////////////////////////////////////////////////////////////////////////////
// Reference;works reasonably well
// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_VE(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega, D)
    const Rcpp::List & configuration // List of config values
) {

  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);   // (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);   // (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);   // (n,p)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (d,p)
  const arma::rowvec & dm = Rcpp::as<arma::rowvec>(params["dm"]); // (p,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const auto M0    = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const auto S0    = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const auto Mu0   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
  const auto Delta0= Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
  const auto Tau0  = Rcpp::as<arma::mat>(params["Tau"]); // (q,p)
  const auto metadata = tuple_metadata(M0, S0, Mu0, Delta0);
  enum { M_ID, S_ID, Mu_ID, Delta_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<M_ID>(parameters.data()) = M0;
  metadata.map<S_ID>(parameters.data()) = S0;
  metadata.map<Mu_ID>(parameters.data()) = Mu0;
  metadata.map<Delta_ID>(parameters.data()) = Delta0;
  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());

  arma::mat Tau = Tau0 ;
  const arma::mat XB = X * B ;
  const arma::rowvec domega = diagvec(Omega).t() ;

  // Optimize
  auto objective_and_grad = [&metadata, &Y, &X, &O, &Tau, &Omega, &domega, &dm, &XB](const double * params, double * grad) -> double {
    const double n = double(X.n_rows) ;
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    const arma::mat Mu = metadata.map<Mu_ID>(params);
    const arma::mat Delta = metadata.map<Delta_ID>(params);

    arma::mat Delta2 = Delta % Delta ;
    arma::mat MumXB = Mu - XB ;
    arma::mat S2 = S % S ;
    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
    arma::mat A2 = trunc_exp(M + .5 * S2) ;
    arma::mat A  = A1 % (A2 * Tau) ;
    arma::mat A_T = A2 % (A1 * Tau.t()) ;
    arma::mat Sigma = (M.t() * M + diagmat(sum(S2,0))) / n ;
    arma::rowvec d = mean(pow(MumXB, 2) + Delta2, 0) ;

    double objective =
      accu(A - Y % (O + Mu + M * Tau))
      - .5 * (accu(log(Delta2)) + accu(log(S2)))
      + .5 * n * accu(log(d)) // profiled in dm
      + .5 * n * (accu(Omega % Sigma) + log_det_sympd(Sigma)) ;

    metadata.map<Mu_ID>(grad)    = (MumXB.each_row() / d) + A - Y ;
    metadata.map<Delta_ID>(grad) = (Delta.each_row() / d) + Delta % A - pow(Delta, -1) ;
    metadata.map<M_ID>(grad)     = M * Omega + A_T - Y * Tau.t() ;
    metadata.map<S_ID>(grad)     = (S.each_row() % domega) + S % A_T - pow(S, -1) ;

    return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

  arma::mat M = metadata.copy<M_ID>(parameters.data());
  arma::mat S = metadata.copy<S_ID>(parameters.data());
  arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
  arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
  arma::mat Delta2 = Delta % Delta ;
  arma::mat MumXB = Mu - XB ;
  arma::mat S2 = S % S ;
  arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
  arma::mat A2 = trunc_exp(M + .5 * S2) ;
  arma::mat A  = A1 % (A2 * Tau) ;

  arma::vec log_pi = trunc_log(mean(Tau,1)) ;
  if (Tau.n_rows > 1) {
    Tau = (M.t() * Y - A2.t() * A1) ;
    Tau.each_col() += log_pi - 1 ;
    Tau.each_col( [](arma::vec& x){
      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
    }) ;
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = static_cast<int>(result.status),
    Rcpp::Named("iterations") = result.nb_iterations,
    Rcpp::Named("objective") = result.objective,
    Rcpp::Named("M") = M,
    Rcpp::Named("S") = S,
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Delta") = Delta,
    Rcpp::Named("Tau") = Tau,
    Rcpp::Named("A") = A
  );
}

// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_VE_blocks(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//   const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
//   const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const auto M0   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const auto S0   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const auto Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const auto Delta   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const auto metadata = tuple_metadata(M0, S0);
//   enum { M_ID, S_ID }; // Names for metadata indexes
//
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<M_ID>(parameters.data()) = M0;
//   metadata.map<S_ID>(parameters.data()) = S0;
//
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   const arma::mat XB = X * B ;
//   const arma::rowvec domega = diagvec(Omega).t() ;
//   const arma::mat YT = Y * T.t() ;
//   const arma::mat A1 = trunc_exp(O + Mu + .5 * Delta % Delta) ;
//
//   // Optimize
//   auto objective_and_grad = [&metadata, &Y, &YT, &X, &O, &T, &Mu, &Delta, &A1, &Omega, &domega, &XB](const double * params, double * grad) -> double {
//     const double n = double(X.n_rows) ;
//     const arma::mat M = metadata.map<M_ID>(params);
//     const arma::mat S = metadata.map<S_ID>(params);
//
//     arma::mat S2 = S % S ;
//     arma::mat A2 = trunc_exp(M + .5 * S2) ;
//     arma::mat A  = A1 % (A2 * T) ;
//     arma::mat A_T = A2 % (A1 * T.t()) ;
//     arma::mat Sigma = (M.t() * M + diagmat(sum(S2,0))) / n ;
//
//     double objective =
//       accu(A - Y % (O + Mu + M * T))
//       + .5 * n * (
//           accu(Omega % Sigma) + log_det_sympd(Sigma)
//       ) - 0.5 * accu(log(S2)) ;
//
//     metadata.map<S_ID>(grad)     = (S.each_row() % domega) + S % A_T - pow(S, -1) ;
//     metadata.map<M_ID>(grad)     = M * Omega + A_T - YT ;
//
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   arma::mat M = metadata.copy<M_ID>(parameters.data());
//   arma::mat S = metadata.copy<S_ID>(parameters.data());
//
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("M") = M,
//     Rcpp::Named("S") = S);
// }
//
// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_VE_species(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//   //std::cout << "optim_plnblockbis_VE" << std::endl;
//   const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
//   const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
//   const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const auto M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const auto S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const auto Mu0   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const auto Delta0   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const auto metadata = tuple_metadata(Mu0, Delta0);
//   enum { Mu_ID, Delta_ID }; // Names for metadata indexes
//
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<Mu_ID>(parameters.data()) = Mu0;
//   metadata.map<Delta_ID>(parameters.data()) = Delta0;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   const arma::mat XB = X * B ;
//   const arma::rowvec dm = 1. / diagvec(D).t() ;
//   const arma::mat A2 = trunc_exp(M + .5 * (S % S)) ;
//
//   // Optimize
//   auto objective_and_grad = [&metadata, &Y, &X, &O, &T, &Omega, &M, &S, &dm, &XB, &A2](const double * params, double * grad) -> double {
//     const double n = double(X.n_rows) ;
//     const arma::mat Mu = metadata.map<Mu_ID>(params);
//     const arma::mat Delta = metadata.map<Delta_ID>(params);
//
//     arma::mat Delta2 = Delta % Delta ;
//     arma::mat MumXB = Mu - XB ;
//     arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
//     arma::mat A  = A1 % (A2 * T) ;
//     arma::mat A_T = A2 % (A1 * T.t()) ;
//     arma::rowvec d = mean(MumXB % MumXB + Delta2, 0) ;
//
//     double objective =
//       accu(A - Y % (O + Mu + M * T))
//       + .5 * n * (accu(dm % d) + accu(log(d)))
//       - 0.5 * accu(log(Delta2)) ;
//
//     metadata.map<Delta_ID>(grad) = (Delta.each_row() / d) + Delta % A - pow(Delta, -1) ;
//     metadata.map<Mu_ID>(grad)    = (MumXB.each_row() / d) + A - Y ;
//
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
//   arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
//
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("Mu") = Mu,
//     Rcpp::Named("Delta") = Delta);
// }

// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_M(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//
//   // OPTIMIZER DECLARATION
//   const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const auto metadata = tuple_metadata(init_M);
//   enum { M_ID, }; // Names for metadata indexes
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<M_ID>(parameters.data()) = init_M;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   // QUANTITY THAT REMAIN FIXED DURING OPTIMIZATION
//   const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const arma::mat & S     = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const arma::mat & Mu    = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const arma::mat A1  = trunc_exp(O + Mu + .5 * Delta % Delta) ;
//   const arma::mat A1T = A1 * T.t() ;
//   const arma::mat YT  = Y * T.t() ;
//   const arma::mat S2  = S % S ;
//
//   // FUNCTION THAT COMPUTE OBJECTIVE AND GRADIENT
//   auto objective_and_grad = [&metadata, &Y, &T, &YT, &S2, &A1, &A1T, &Omega](const double * params, double * grad) -> double {
//     const arma::mat M = metadata.map<M_ID>(params);
//     arma::mat A2 = trunc_exp(M + .5 * S2) ;
//     arma::mat A  = A1 % (A2 * T) ;
//     arma::mat AT = A2 % A1T ;
//
//     double objective = accu(A - Y % (M * T)) + .5 * accu((M * Omega) % M) ;
//     metadata.map<M_ID>(grad) = M * Omega + AT - YT  ;
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   // SENDING RESULTS RESULTS
//   arma::mat M = metadata.copy<M_ID>(parameters.data());
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("M") = M);
// }
//
// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_S(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//
//   // OPTIMIZER DECLARATION
//   const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,p)
//   const auto metadata = tuple_metadata(init_S);
//   enum { S_ID }; // Names for metadata indexes
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<S_ID>(parameters.data()) = init_S;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   // QUANTITY THAT REMAIN FIXED DURING OPTIMIZATION
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const arma::mat & M     = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const arma::mat & Mu    = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const arma::mat & T     = Rcpp::as<arma::mat>(params["T"]); // (n,q)
//   const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const arma::mat A1 = trunc_exp(O + Mu + .5 * Delta % Delta) ;
//   const arma::mat A1T = A1 * T.t() ;
//   const arma::vec d_Omega = diagvec(Omega) ;
//
//   // FUNCTION THAT COMPUTE OBJECTIVE AND GRADIENT
//   auto objective_and_grad = [&metadata, &T, &d_Omega, &M, &A1, &A1T](const double * params, double * grad) -> double {
//     const arma::mat S = metadata.map<S_ID>(params);
//     arma::mat S2 = S % S ;
//     arma::mat A2 = trunc_exp(M + .5 * S2) ;
//
//     double objective = accu(A1 % (A2 * T)) + .5 * accu(S2 * d_Omega) - .5 * accu(log(S2)) ;
//     metadata.map<S_ID>(grad) = S % (A1T % A2) + S.each_row() % d_Omega.t() - 1./S  ;
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   // SENDING RESULTS
//   arma::mat S = metadata.copy<S_ID>(parameters.data());
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("S") = S);
// }
//
// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_Mu(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//
//   // OPTIMIZER DECLARATION
//   const auto init_Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
//   const auto metadata = tuple_metadata(init_Mu);
//   enum { Mu_ID }; // Names for metadata indexes
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<Mu_ID>(parameters.data()) = init_Mu;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   // QUANTITY THAT REMAIN FIXED DURING OPTIMIZATION
//   const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // responses (n,d)
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
//   const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (d,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const arma::mat & M     = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const arma::mat & S     = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const arma::mat & T     = Rcpp::as<arma::mat>(params["T"]); // (q,p)
//   const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const arma::mat XB  = X * B ;
//   const arma::vec dm  = 1./arma::diagvec(D) ;
//   const arma::mat expDA2T  = trunc_exp(O + .5 * Delta % Delta) % (trunc_exp(M + .5 * S % S) * T) ;
//
//   // FUNCTION THAT COMPUTE OBJECTIVE AND GRADIENT
//   auto objective_and_grad = [&metadata, &Y, &dm, &XB, &expDA2T](const double * params, double * grad) -> double {
//     const arma::mat Mu = metadata.map<Mu_ID>(params);
//     arma::mat A  = trunc_exp(Mu) % expDA2T ;
//     arma::mat MumXB = Mu - XB ;
//
//     double objective = accu(A - Y % Mu) + .5 * accu(pow(MumXB, 2)*dm) ;
//     metadata.map<Mu_ID>(grad) = MumXB.each_row() % dm.t() + A - Y  ;
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   // SENDING RESULTS RESULTS
//   arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("Mu") = Mu);
// }
//
// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_Delta(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//
//   // OPTIMIZER DECLARATION
//   const auto init_Delta   = Rcpp::as<arma::mat>(params["Delta"]);
//   const auto metadata = tuple_metadata(init_Delta);
//   enum { Delta_ID }; // Names for metadata indexes
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<Delta_ID>(parameters.data()) = init_Delta;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   // QUANTITY THAT REMAIN FIXED DURING OPTIMIZATION
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (n,q)
//   const arma::mat & Mu= Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (n,q)
//   const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const arma::mat expMuA2T  = trunc_exp(O + Mu) % (trunc_exp(M + .5 * S % S) * T)  ;
//   const arma::vec dm  = 1./arma::diagvec(D) ;
//
//   // FUNCTION THAT COMPUTE OBJECTIVE AND GRADIENT
//   auto objective_and_grad = [&metadata, &T, &dm, &expMuA2T](const double * params, double * grad) -> double {
//     const arma::mat Delta = metadata.map<Delta_ID>(params);
//     arma::mat Delta2 = Delta % Delta ;
//     arma::mat A  = exp(.5 * Delta2) % expMuA2T   ;
//
//     double objective = accu(A) + .5 * accu(Delta2 * dm) - .5 * accu(log(Delta2)) ;
//     metadata.map<Delta_ID>(grad) = Delta % A  + Delta.each_row() % dm.t() - 1./Delta  ;
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   // SENDING RESULTS
//   arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("Delta") = Delta);
// }


// // [[Rcpp::export]]
// Rcpp::List optim_plnblockbis_VE_profiled(
//     const Rcpp::List & data  , // List(Y, X, O, w)
//     const Rcpp::List & params, // List(M, S, T, B, Omega, D)
//     const Rcpp::List & configuration // List of config values
// ) {
//   const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
//   const arma::mat & XtXm = Rcpp::as<arma::mat>(data["XtXm"]); // covariates (n,d)
//   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
//   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
//   const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
//   const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
//   const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
//   const auto init_M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//   const auto init_S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
//   const auto init_Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
//   const auto init_Delta   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
//   const auto metadata = tuple_metadata(init_M, init_S, init_Mu, init_Delta);
//   enum { M_ID, S_ID, Mu_ID, Delta_ID }; // Names for metadata indexes
//
//   auto parameters = std::vector<double>(metadata.packed_size);
//   metadata.map<M_ID>(parameters.data()) = init_M;
//   metadata.map<S_ID>(parameters.data()) = init_S;
//   metadata.map<Mu_ID>(parameters.data()) = init_Mu;
//   metadata.map<Delta_ID>(parameters.data()) = init_Delta;
//   auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
//
//   double w_bar = accu(w);
//   const arma::mat ProjOrthX = arma::eye(X.n_rows, X.n_rows) - X * XtXm * X.t() ;
//
//   // Optimize
//   auto objective_and_grad = [&metadata, &Y, &ProjOrthX, &O, &w, &w_bar, &T](const double * params, double * grad) -> double {
//     const arma::mat M = metadata.map<M_ID>(params);
//     const arma::mat S = metadata.map<S_ID>(params);
//     const arma::mat Mu = metadata.map<Mu_ID>(params);
//     const arma::mat Delta = metadata.map<Delta_ID>(params);
//
//     arma::mat MumXB = ProjOrthX * Mu  ;
//     arma::mat Delta2 = Delta % Delta ;
//     arma::mat S2 = S % S ;
//     arma::rowvec dm = pow(w.t() * (pow(MumXB, 2) + Delta2) / w_bar, -1);
//     arma::mat Omega = w_bar * inv_sympd(M.t() * (M.each_col() % w) + diagmat(w.t() * S2));
//     arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
//     arma::mat A2 = trunc_exp(M + .5 * S2) ;
//     arma::mat A  = A1 % (A2 * T) ;
//     arma::mat A_T = A2 % (A1 * T.t()) ;
//
//     double objective = accu(w.t() * (A - Y % (O + Mu + M * T)))
//       - 0.5 * (accu(w.t() * log(S2)) + accu(w.t() * log(Delta2)))
//       - 0.5 * w_bar * (log_det_sympd(Omega) + accu(log(dm)))
//       ;
//
//     metadata.map<S_ID>(grad)     = diagmat(w) *  (S.each_row() % diagvec(Omega).t() + S % A_T - pow(S, -1)) ;
//     metadata.map<Delta_ID>(grad) = diagmat(w) * ((Delta.each_row() % dm) + Delta % A - pow(Delta, -1)) ;
//     metadata.map<M_ID>(grad)  = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;
//     metadata.map<Mu_ID>(grad) = diagmat(w) * ((MumXB.each_row() % dm) + A - Y) ;
//
//     return objective;
//   };
//   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
//
//   arma::mat M = metadata.copy<M_ID>(parameters.data());
//   arma::mat S = metadata.copy<S_ID>(parameters.data());
//   arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
//   arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
//
//   return Rcpp::List::create(
//     Rcpp::Named("status") = static_cast<int>(result.status),
//     Rcpp::Named("iterations") = result.nb_iterations,
//     Rcpp::Named("objective") = result.objective,
//     Rcpp::Named("M") = M,
//     Rcpp::Named("S") = S,
//     Rcpp::Named("Mu") = Mu,
//     Rcpp::Named("Delta") = Delta
//   );
// }
