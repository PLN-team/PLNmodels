// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List fn_optim_PLNnetwork_Cpp(const arma::vec par,
                               double log_detOmega,
                               const arma::mat Omega,
                               const arma::mat Y,
                               const arma::mat X,
                               const arma::mat O,
                               double KY) {

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta = par.subvec(0      , p*d      -1) ; Theta.reshape(p,d) ;
  arma::mat M     = par.subvec(p*d    , p*(n+d) - 1) ; M.reshape(n,p) ;
  arma::mat S     = par.subvec(p*(n+d), p*(2*n+d)-1) ; S.reshape(n,p) ;

  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double objective = accu(A - Y % Z - .5*log(S)) -.5*(n*log_detOmega + n*p - trace(Omega*(diagmat(sum(S, 0)) + M.t() * M))) + KY ;

  arma::vec grd_Theta = vectorise(X.t() * (A-Y));
  arma::vec grd_M     = vectorise(M * Omega + A-Y) ;
  arma::vec grd_S     = vectorise(.5 * (ones(n) * diagvec(Omega).t() + A - 1/S));

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_M),grd_S) ;

  return Rcpp::List::create(Rcpp::Named("objective") = objective,
                            Rcpp::Named("gradient" ) = grad);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
load("~/svn/sparsepca/Data/OTU/oaks.RData")

formula <- Data$count[1:10, 1:5] ~ 1 + offset(log(Data$offset[1:10, 1:5]))

frame  <- model.frame(formula)
Y      <- model.response(frame)
X      <- model.matrix(formula)
O      <- model.offset(frame)
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
KY <- 1 ## constant quantity in the objective

## declare the objective and gradient functions for optimization
fn_optim_PLNnetwork <- function(par,logDetOmega,Omega,Y,X,O) {

      Theta <- matrix(par[1:(p*d)]                         , p,d)
      M     <- matrix(par[p*d          + 1:(n*p)], n,p)
      S     <- matrix(par[(n+d)*p + 1:(n*p)], n,p)

      Z <- O + tcrossprod(X, Theta) + M
      A <- exp (Z + .5*S)

      logP.Z  <- n/2 * (logDetOmega - sum(diag(Omega)*colMeans(S))) - .5*sum(diag(Omega %*% crossprod(M)))

      gr.Theta <- crossprod(X, A - Y)
      gr.M  <- M %*% Omega + A - Y
      gr.S  <- .5 * (matrix(rep(diag(Omega),n), n, p, byrow = TRUE) + A - 1/S)

      return(list(
        "objective" = sum(as.numeric(A - Y*Z)) - logP.Z - .5*sum(log(S)+1) + KY,
        "gradient"  = c(gr.Theta,gr.M,gr.S)
      ))
  }

logLM  <- lapply(1:ncol(Y), function(j) lm.fit(X, log(1 + Y[, j]), offset = O[,j]))
Theta <- do.call(rbind, lapply(logLM, coefficients))
M <- sapply(logLM, residuals)
S <- matrix(1e-3, n, p)

Omega <- solve(cov(M))
logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus

par0 <- c(Theta, M, S)

outR <- fn_optim_PLNnetwork(par0, logDetOmega, Omega, Y, X, O)
outC <- fn_optim_PLNnetwork_Cpp(par0, logDetOmega, Omega, Y, X, O, KY)

library(microbenchmark)
res <- microbenchmark(outR = fn_optim_PLNnetwork(par0, logDetOmega, Omega, Y, X, O),
                      outC = fn_optim_PLNnetwork_Cpp(par0, logDetOmega, Omega, Y, X, O, KY))
print(summary(res))
ggplot2::autoplot(res)

## checking derivative
library(nloptr)
# out <- check.derivatives(par0, func = function(x) {
#                 fn_optim_PLNnetwork(x, logDetOmega, Omega, Y, X, O)$objective
#             }, func_grad = function(x) {
#                 fn_optim_PLNnetwork(x, logDetOmega, Omega, Y, X, O)$gradient
#             })
out <- check.derivatives(par0, func = function(x) {
                fn_optim_PLNnetwork_Cpp(x, logDetOmega, Omega, Y, X, O, KY)$objective
            }, func_grad = function(x) {
                fn_optim_PLNnetwork_Cpp(x, logDetOmega, Omega, Y, X, O, KY)$gradient
            }, check_derivatives_print = "all")

*/
