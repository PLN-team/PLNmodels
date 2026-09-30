#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "graphical_lasso.h"

// Backend of the R function graphical_lasso(), which validates the input and
// expands a scalar penalty to a matrix before calling this. The stall_patience
// default here must track graphical_lasso::kStallPatience (Rcpp attributes cannot
// parse a namespaced constant as a default).
// [[Rcpp::export]]
Rcpp::List cpp_graphical_lasso(const arma::mat & S, const arma::mat & rho,
                               double thr, int maxit,
                               Rcpp::Nullable<Rcpp::NumericMatrix> w_init = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> wi_init = R_NilValue,
                               bool trace = false,
                               int stall_patience = 1000) {
  graphical_lasso::State warm;
  if (w_init.isNotNull() && wi_init.isNotNull()) {
    warm.W = Rcpp::as<arma::mat>(w_init.get());
    warm.X = Rcpp::as<arma::mat>(wi_init.get());
    warm.filled = true;
  }
  graphical_lasso::Result res = graphical_lasso::solve(S, rho, thr, maxit, warm.filled ? &warm : nullptr, trace, stall_patience);
  return Rcpp::List::create(
    Rcpp::Named("w")         = res.W,
    Rcpp::Named("wi")        = res.X,
    Rcpp::Named("niter")     = res.niter,
    Rcpp::Named("converged") = res.converged,
    Rcpp::Named("status")    = graphical_lasso::status_name(res.status),
    Rcpp::Named("delta")     = res.delta,
    Rcpp::Named("shift")     = res.shift,
    Rcpp::Named("dw_trace")  = res.dw_trace
  );
}
