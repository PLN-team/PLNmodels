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
  //std::cout << "plnblockbis_vloglik" << std::endl;
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  //const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat &T  = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword p = Y.n_cols;
  const arma::uword q = M.n_cols;
  const arma::uword n = M.n_rows;
  const arma::mat S2 = S % S ;
  const arma::mat Mu_XB = Mu - X * B ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat Z = O + Mu + M * T;
  const arma::mat A = trunc_exp(O + Mu + .5 * Delta2) % (trunc_exp(M + .5 * S2) * T) ;
  const arma::vec Dm = 1. / arma::diagvec(D) ;

  // Element-wise log-likelihood
  return(
         sum(Y % Z - A, 1) - logfact(Y)                                      // Poisson term
      + 0.5 * real(log_det(Omega)) + 0.5 * sum(log(S2), 1) + 0.5 * double(q) // Gaussian-block term
      - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)                 // ...
      - 0.5 * real(log_det(D)) + 0.5 * sum(log(Delta2), 1) + 0.5 * double(p) // Gaussian-species term
      - 0.5 * (pow(Mu_XB, 2) + Delta2) * Dm                                  // ...
      + accu(log_pi.t() * T) - accu(T % arma::trunc_log(T))                  // Multinomial-block term
  ) ;
}

// [[Rcpp::export]]
double plnblockbis_loglik(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  //std::cout << "plnblockbis_loglik" << std::endl;
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (n,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (p,p)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));

  const arma::uword p = Y.n_cols;
  const arma::uword q = M.n_cols;
  const arma::uword n = M.n_rows;
  const arma::mat S2 = S % S ;
  const arma::mat MumXB = Mu - X * B ;
  const arma::mat Delta2 = Delta % Delta ;
  const arma::mat Z = O + Mu + M * T;
  const arma::mat A = trunc_exp(O + Mu + .5 * Delta2) % (trunc_exp(M + .5 * S2) * T) ;
  const arma::vec Dm = 1. / arma::diagvec(D) ;

  // Element-wise log-likelihood
  return (as_scalar(w.t() *
          (sum(Y % Z - A, 1) - logfact(Y)
            + 0.5 * real(log_det(Omega)) + 0.5 * sum(log(S2), 1) + 0.5 * double(q)
            - 0.5 * sum( (M * Omega) % M + S2 * diagmat(Omega), 1)
            - 0.5 * real(log_det(D)) + 0.5 * sum(log(Delta2), 1) + 0.5 * double(p)
            - 0.5 * (pow(MumXB, 2) + Delta2) * Dm
          )) + accu(log_pi.t() * T) - accu(T % arma::trunc_log(T))
    );
}

// [[Rcpp::export]]
Rcpp::List  optim_plnblockbis_Omega(
  const arma::mat & M, // (n,q)
  const arma::mat & S, // (n,q)
  const arma::vec & w  // (n)
) {
  double w_bar = accu(w);
  arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * (S % S)))  ;

  //std::cout << "inv_sympd 2" << std::endl;
  const arma::uword n = M.n_rows;
  return Rcpp::List::create(
    Rcpp::Named("Sigma") = (1./w_bar) * nSigma,
    Rcpp::Named("Omega") = w_bar * inv_sympd(nSigma));
}


// [[Rcpp::export]]
Rcpp::List  optim_plnblockbis_B(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,  // List(M, S, T, B, Omega)
    const Rcpp::List & configuration // (n)
) {
  //std::cout << "optim_plnblockbis_B" << std::endl;
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  //double w_bar = accu(w);

  //std::cout<< "inv_sympd 3" << std::endl;
  arma::mat B = inv_sympd((X.t() * diagmat(w) * X)) * X.t()  * diagmat(w) * Mu;
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
  const arma::mat MumXB = Mu - X * B ;
  arma::rowvec d = w.t() * (MumXB % MumXB + Delta % Delta) / accu(w);

  return Rcpp::List::create(
    Rcpp::Named("D") = diagmat(d.t())
  ) ;
}

// [[Rcpp::export]]
arma::mat optim_plnblockbis_Tau(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params  // List(M, S, T, B, Omega)
) {
  //std::cout << "optim_plnblockbis_Tau" << std::endl;
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const arma::vec log_pi = arma::trunc_log(mean(T,1));
  const arma::mat & Mu = Rcpp::as<arma::mat>(params["Mu"]); // (n,p)
  const arma::mat & Delta = Rcpp::as<arma::mat>(params["Delta"]); // (n,p)
  const arma::mat A1 = trunc_exp(O + Mu + .5* Delta % Delta) ;
  const arma::mat A2 = trunc_exp(M + .5 * S % S ) ;
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)

  arma::mat Tau = M.t() * diagmat(w) * Y - A2.t()* diagmat(w) * A1  ;
  Tau.each_col() += log_pi ;
  Tau.each_col( [](arma::vec& x){
    x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
  }) ;
  return Tau ;
}


// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_VE_blocks(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega, D)
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
  const auto Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
  const auto Delta   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
  const auto metadata = tuple_metadata(init_M, init_S);
  enum { M_ID, S_ID}; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<M_ID>(parameters.data()) = init_M;
  metadata.map<S_ID>(parameters.data()) = init_S;

  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
  set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));

  const arma::mat XB = X * B ;

  // Optimize
  auto objective_and_grad = [&metadata, &Mu, &Delta, &Y, &X, &O, &T, &Omega, &w, &XB](const double * params, double * grad) -> double {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);

    double w_bar = accu(w);
    arma::mat Delta2 = Delta % Delta ;
    arma::mat MumXB = Mu - XB ;
    arma::mat S2 = S % S ;
    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
    arma::mat A2 = trunc_exp(M + .5 * S2) ;
    arma::mat A  = A1 % (A2 * T) ;
    arma::mat A_T = A2 % (A1 * T.t()) ;
    arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) ;
    arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;
    arma::rowvec omega = diagvec(Omega).t() ;

    double objective =
      accu(w.t() * (A - Y % (O + Mu + M * T)))  // ELBO term for Poisson
      - 0.5 * accu(w.t() * log(S2))             // ELBO term for blocks
      + .5 * trace(Omega * nSigma) ;            // ...

      metadata.map<S_ID>(grad) = diagmat(w) *  (S.each_row() % omega + S % A_T - pow(S, -1)) ; // question omega ici ?
      metadata.map<M_ID>(grad) = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;

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


// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_VE_species(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega, D)
    const Rcpp::List & configuration // List of config values
) {
  //std::cout << "optim_plnblockbis_VE" << std::endl;
  const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
  const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
  const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
  const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (n,q)
  const arma::mat & T = Rcpp::as<arma::mat>(params["T"]); // (q,p)
  const arma::mat & Omega = Rcpp::as<arma::mat>(params["Omega"]); // (q,q)
  const auto M   = Rcpp::as<arma::mat>(params["M"]); // (n,q)
  const auto S   = Rcpp::as<arma::mat>(params["S"]); // (n,q)
  const auto init_Mu   = Rcpp::as<arma::mat>(params["Mu"]); // (n,q)
  const auto init_Delta   = Rcpp::as<arma::mat>(params["Delta"]); // (n,q)
  const auto metadata = tuple_metadata(init_Mu, init_Delta);
  enum { Mu_ID, Delta_ID }; // Names for metadata indexes

  auto parameters = std::vector<double>(metadata.packed_size);
  metadata.map<Mu_ID>(parameters.data()) = init_Mu;
  metadata.map<Delta_ID>(parameters.data()) = init_Delta;
  auto optimizer = new_nlopt_optimizer(configuration, parameters.size());
  set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));

  const arma::mat XB = X * B ;
  double w_bar = accu(w);

  // Optimize
  auto objective_and_grad = [&metadata, &M, &S, &Y, &X, &O, &T, &Omega, &w, &w_bar, &XB](const double * params, double * grad) -> double {
    const arma::mat Mu = metadata.map<Mu_ID>(params);
    const arma::mat Delta = metadata.map<Delta_ID>(params);

    arma::mat Delta2 = Delta % Delta ;
    arma::mat MumXB = Mu - XB ;
    arma::mat S2 = S % S ;
    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
    arma::mat A2 = trunc_exp(M + .5 * S2) ;
    arma::mat A  = A1 % (A2 * T) ;
    arma::mat A_T = A2 % (A1 * T.t()) ;
    arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;
    arma::rowvec omega = 1./diagvec(Omega).t() ;

    double objective =
      accu(w.t() * (A - Y % (O + Mu + M * T)))  // ELBO term for Poisson
      - 0.5 * accu(w.t() * log(Delta2))         // ELBO term for species
      + 0.5 * w_bar * accu(log(d)) ;            // ...

      metadata.map<Delta_ID>(grad) = diagmat(w) * ((Delta.each_row() / d) + Delta % A - pow(Delta, -1)) ;
      metadata.map<Mu_ID>(grad) = diagmat(w) * ((MumXB.each_row() / d) + A - Y) ;

      return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

  arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
  arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
  arma::rowvec d = w.t() * (pow(Mu - XB,2) + Delta % Delta) / w_bar;

  return Rcpp::List::create(
    Rcpp::Named("status") = static_cast<int>(result.status),
    Rcpp::Named("iterations") = result.nb_iterations,
    Rcpp::Named("objective") = result.objective,
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Delta") = Delta,
    Rcpp::Named("D") = arma::diagmat(d));
}






////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List optim_plnblockbis_VE(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(M, S, T, B, Omega, D)
    const Rcpp::List & configuration // List of config values
) {
  //std::cout << "optim_plnblockbis_VE" << std::endl;
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

  const arma::mat XB = X * B ;
  double w_bar = accu(w);
  arma::vec log_pi = arma::trunc_log(mean(T,1));

  // Optimize
  auto objective_and_grad = [&metadata, &Y, &X, &O, &T, &Omega, &w, &w_bar, &XB, &log_pi](const double * params, double * grad) -> double {
    const arma::mat M = metadata.map<M_ID>(params);
    const arma::mat S = metadata.map<S_ID>(params);
    const arma::mat Mu = metadata.map<Mu_ID>(params);
    const arma::mat Delta = metadata.map<Delta_ID>(params);

    arma::mat Delta2 = Delta % Delta ;
    arma::mat MumXB = Mu - XB ;
    arma::mat S2 = S % S ;
    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;
    arma::mat A2 = trunc_exp(M + .5 * S2) ;

    // Changement ici : ajout du facteur  (1/w_bar)
    //arma::mat Tau = (1/w_bar) * (M.t() * diagmat(w) * Y - A2.t()* diagmat(w) * A1)  ;
    arma::mat Tau = (M.t() * diagmat(w) * Y - A2.t()* diagmat(w) * A1)  ;
    Tau.each_col() += log_pi ;
    Tau.each_col( [](arma::vec& x){
      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
    }) ;

    arma::mat A  = A1 % (A2 * T) ;
    arma::mat A_T = A2 % (A1 * T.t()) ;
    arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) ;
    arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;
    arma::rowvec omega = diagvec(Omega).t() ;

    double objective =
      accu(w.t() * (A - Y % (O + Mu + M * T)))  // ELBO term for Poisson
      - 0.5 * accu(w.t() * log(S2))             // ELBO term for blocks
      + .5 * trace(Omega * nSigma)              // ...
      - 0.5 * accu(w.t() * log(Delta2))         // ELBO term for species
      + 0.5 * w_bar * accu(log(d)) ;            // ...

    metadata.map<S_ID>(grad)     = diagmat(w) *  ((S.each_row() % omega) + S % A_T - pow(S, -1)) ; // question omega ici ?
    metadata.map<Delta_ID>(grad) = diagmat(w) * ((Delta.each_row() / d) + Delta % A - pow(Delta, -1)) ;
    metadata.map<M_ID>(grad)  = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;
    metadata.map<Mu_ID>(grad) = diagmat(w) * ((MumXB.each_row() / d) + A - Y) ;

    return objective;
  };
  OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());
    arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());
    arma::rowvec d = w.t() * (pow(Mu - XB,2) + Delta % Delta) / w_bar;
    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta%Delta) ;
    arma::mat A2 = trunc_exp(M + .5 * S%S) ;
    arma::mat Tau = M.t() * diagmat(w) * Y - A2.t() * diagmat(w) * A1  ;
    Tau.each_col() += log_pi ;
    Tau.each_col( [](arma::vec& x){
      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;
    }) ;
    arma::mat A  = A1 % (A2 * T) ;

    return Rcpp::List::create(
      Rcpp::Named("status") = static_cast<int>(result.status),
      Rcpp::Named("iterations") = result.nb_iterations,
      Rcpp::Named("objective") = result.objective,
      Rcpp::Named("M") = M,
      Rcpp::Named("S") = S,
      Rcpp::Named("Mu") = Mu,
      Rcpp::Named("Delta") = Delta,
      Rcpp::Named("Tau") = Tau,
      Rcpp::Named("D") = diagmat(d),
      Rcpp::Named("A") = A
    );
}
