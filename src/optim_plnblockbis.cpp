// BESOIN DE REVERIFIER TOUTES LES EXPRESSIONS D'ELBO ET ASSIMIL´E
// AINSI QUE LES EXPRESSIONS DE GRADIENTS
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
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat &T  = Rcpp::as<arma::mat>(params["T"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword q = M.n_cols;
  const arma::mat S2 = S % S ;
  const arma::mat XB = X * B ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat mr = O + Mu ;
  const arma::mat Z = mr + M * T;
  const arma::mat A = trunc_exp(mr + .5 * Delta) % (trunc_exp(M + .5 * S2) * T) ;
  const arma::mat Dm = inv_sympd(D);

  // Element-wise log-likelihood
  // rq: j'ai laissé tomber les constantes
  return(
    0.5 * real(log_det(Omega)) - logfact(Y) + sum(Y % Z - A, 1)
    - 0.5 * real(log_det(D)) - accu(T % arma::trunc_log(T))
    + accu(log_pi.t() * T)
    + 0.5 * sum(log(S2), 1) + 0.5 * sum(log(Delta2), 1)
    - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)
    - 0.5 * sum( (Mu * Omega) % Mu + Delta2 * Dm, 1)
    - 0.5 * sum(XB % Dm % XB.t(), 1)
    + sum(XB % Dm % (XB.t()), 1)
    + sum(Mu % Dm % (XB.t()), 1)
  ) ;

}

// [[Rcpp::export]]
arma::vec plnblockbis_loglik(
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
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat &T  = Rcpp::as<arma::mat>(params["T"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword q = M.n_cols;
  const arma::mat S2 = S % S ;
  const arma::mat XB = X * B ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat mr = O + Mu ;
  const arma::mat Z = mr + M * T;
  const arma::mat A = trunc_exp(mr + .5 * Delta2) % (trunc_exp(M + .5 * S2) * T) ;
  const arma::mat Dm = inv_sympd(D);

  // Element-wise log-likelihood
  // rq: j'ai laissé tomber les constantes
  return(
    0.5 * real(log_det(Omega)) - logfact(Y) + sum(Y % Z - A, 1)
    - 0.5 * real(log_det(D)) - accu(T % arma::trunc_log(T))
    + accu(log_pi.t() * T)
    + 0.5 * sum(log(S2), 1) + 0.5 * sum(log(Delta2), 1)
    - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)
    - 0.5 * sum( (Mu * Omega) % Mu + Delta2 * Dm, 1)
    - 0.5 * sum(XB % Dm % XB.t(), 1)
    + sum(XB % Dm % (XB.t()), 1)
    + sum(Mu % Dm % (XB.t()), 1)
  ) ;
}

// [[Rcpp::export]]
Rcpp::List  optim_plnblockbis_Omega(
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
Rcpp::List  optim_plnblockbis_B(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,  // List(M, S, T, B, Omega)
    const arma::vec & w  // (n)
) {
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  double w_bar = accu(w);
  arma::mat B = inv_sympd((X.t() * w * X)) * X.t() * w * Mu;
  return Rcpp::List::create(
    Rcpp::Named("B") = B);
}

// [[Rcpp::export]]
Rcpp::List  optim_plnblockbis_D(
    const arma::mat & X,
    const arma::mat & B,
    const arma::mat & Mu,
    const arma::mat & Delta,
    const arma::vec & w  // (n)
) {
  const arma::mat Mu2 = Mu % Mu ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat XB = X * B ;

  double w_bar = accu(w);
  arma::vec d = sum(Mu2 + Delta2 + XB - 2 * (Mu % XB) , 1);
  arma::mat nD = diagmat(d);

  return Rcpp::List::create(
    Rcpp::Named("D") = (1./w_bar) * nD);
}

// [[Rcpp::export]]
arma::mat optim_plnblockbis_Tau(
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
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat Delta2 = Delta % Delta ;
  const arma::vec log_pi = arma::trunc_log(mean(T,1));
  const arma::mat mr = O + Mu ;
  const arma::mat A1 = trunc_exp(mr + .5 * Delta2) ;
  const arma::mat A2 = trunc_exp(M + .5 * S % S ) ;

  arma::mat Tau = M.t() * Y - A1.t() * A2  ;
  Tau.each_col() += log_pi ;
  Tau.each_col( [](arma::vec& x){
    x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
  }) ;
  return Tau ;
}

// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_VE(
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
  const auto init_Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
  const auto init_Delta   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
  const auto metadata = tuple_metadata(init_M, init_S, init_Mu, init_Delta);
  enum { M_ID, S_ID, Mu_ID, Delta_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<M_ID>(parameters.data()) = init_M;
  metadata.map<S_ID>(parameters.data()) = init_S;
  metadata.map<Mu_ID>(parameters.data()) = init_Mu;
  metadata.map<Delta_ID>(parameters.data()) = init_Delta;
  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
  set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));

  const arma::mat mr = O + X * B ;
  const arma::mat A2 = trunc_exp(mr) ;
  const arma::mat XB = X * B ;

  // Optimize
  auto objective_and_grad = [&metadata, &Y, &X, &O, &mr, &A2,  &T, &Omega, &w, &XB](const double * params, double * grad) -> double {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    const arma::mat Mu = metadata.map<Mu_ID>(params);
    const arma::mat Delta = metadata.map<Delta_ID>(params);

    arma::mat Delta2 = Delta % Delta ;
    arma::mat S2 = S % S ;
    arma::mat Mu2 = Mu % Mu ;
    arma::mat A1  = trunc_exp(O + Mu + .5 * Delta2) ;
    arma::mat A2  = trunc_exp(M + .5 * S2) ;
    arma::mat A   = A1 % (A2 * T) ;
    arma::mat A_T = A2 % (A1 * T.t()) ;
    arma::mat nSigma = M.t() * (M.each_col() % w) + diagmat(w.t() * S2);
    // BESOIN D'INTEGRER w DANS LE CALCUL DE d
    // PEUT-ETRE CONFUSION ENTRE d et d^{-1}
    // SANS DOUTE BESOIN DE GERER LA DIVISION DE D PAR w_bar
    arma::vec d = sum(Mu2 + Delta2 + XB - 2 * (Mu % XB), 0).t();
    arma::mat nD = diagmat(d);
    //
    std::cout << "nD dimensions: " << nD.n_rows << " x " << nD.n_cols << std::endl;
    std::cout << "Mu dimensions: " << Mu.n_rows << " x " << Mu.n_cols << std::endl;
    //
    double w_bar = accu(w);
    arma::mat nDm = diagmat(1./d);
    arma::mat Dm = diagmat((1./w_bar) * 1./d);

    // verifier expression de la fonction objective
    double objective =
      accu(w.t() * (A - Y % (O + Mu + M * T))) - 0.5 * accu(w.t() * log(S2))
      - 0.5 * accu(w.t() * log(Delta2)) + .5 * trace(Omega * nSigma)
      + 0.5 * trace(Mu.t() * (Mu.each_col() % w) * nDm)
      + 0.5 * trace(diagmat(w.t() * Delta2) * nDm);

    metadata.map<M_ID>(grad) = diagmat(w) * (M * Omega + A_T - Y * T.t());
    metadata.map<S_ID>(grad) = diagmat(w) * (S.each_row() % diagvec(Omega).t() + S % A_T - pow(S, -1));
    std::cout << "7!" << std::endl;
    // UNE ERREUR ICI : vérifier l'expression du gradient de Mu, pbm avec A_T
    metadata.map<Mu_ID>(grad) = Y - Mu * nD ;//- A_T
    std::cout << "8!" << std::endl;
    metadata.map<Delta_ID>(grad) = diagmat(w) * (Delta.each_row() % diagvec(Dm).t() + (A2 * T) % A1 % Delta - pow(Delta, -1));

    return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
    arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
    return Rcpp::List::create(
      Rcpp::Named("status") = static_cast<int>(result.status),
      Rcpp::Named("iterations") = result.nb_iterations,
      Rcpp::Named("objective") = result.objective,
      Rcpp::Named("M") = M,
      Rcpp::Named("S") = S,
      Rcpp::Named("Mu") = Mu,
      Rcpp::Named("Delta") = Delta);
}
