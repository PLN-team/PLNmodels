// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List fn_optim_PLN_Cpp(const arma::vec par,
                        const arma::mat Y,
                        const arma::mat X,
                        const arma::mat O,
                        double KY) {

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;  
  
  arma::mat Theta = par.subvec(0      , p*d      -1) ; Theta.reshape(p,d) ;
  arma::mat M     = par.subvec(p*d    , p*(n+d) - 1) ; M.reshape(n,p) ;
  arma::mat S     = par.subvec(p*(n+d), p*(2*n+d)-1) ; S.reshape(n,p) ;

  arma::mat S_bar(sum(S, 0)) ;
  arma::mat MtM = M.t() * M ;
  arma::mat Omega = n * inv_sympd(MtM + diagmat(S_bar));
  double log_detOmega = real(log_det(Omega)) ;
  
  arma::mat Z = O + X * Theta.t() + M;
  arma::mat A = exp (Z + .5 * S) ;

  double logP_Z = .5 * (n * log_detOmega - dot(diagvec(Omega), S_bar) - trace(Omega * MtM))  ;
  double objective = accu(A - Y % Z - .5 * log(S) - .5) - logP_Z + KY ; 
  
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
source("~/git/PLNmodels/R/utils.R")
load("~/svn/sparsepca/Pgm/PCA/Output/Corinne_data.RData")

formula <- Data$count ~ 1 + offset(log(Data$offset))

frame  <- model.frame(formula)
Y      <- model.response(frame)
X      <- model.matrix(formula)
O      <- model.offset(frame)
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
KY <- sum(.logfactorial(Y)) ## constant quantity in the objective

fn_optim_PLN <- function(par,Y,X,O) {
    Theta <- matrix(par[1:(p*d)]                         , p,d)
    M     <- matrix(par[p*d          + 1:(n*p)], n,p)
    S     <- matrix(par[(n+d)*p + 1:(n*p)], n,p)

    Omega <- chol2inv(chol(crossprod(M)/n + diag(colMeans(S))))
    logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus

    Z <- O + tcrossprod(X, Theta) + M
    A <- exp (.trunc(Z + .5*S))
    logP.Z  <- n/2 * (logDetOmega - sum(diag(Omega)*colMeans(S))) - .5*sum(diag(Omega %*% crossprod(M)))

    gr.Theta <- crossprod(X, A - Y)
    gr.M  <- M %*% Omega + A - Y
    gr.S  <- .5 * (matrix(rep(diag(Omega),n), n, p, byrow = TRUE) + A - 1/S)

    return(list(
      "objective" = sum(as.numeric(A - Y*Z)) - logP.Z - .5*sum(log(S)+1) + KY,
      "gradient"  = c(gr.Theta,gr.M,gr.S)
    ))
  }

glmP  <- lapply(1:ncol(Y), function(j) glm.fit(X, Y[, j], offset = O[,j], family = poisson()))
Theta <- do.call(rbind, lapply(glmP, coefficients))
M <- matrix(0, n, p) ## -tcrossprod(covariates,self$init.par$Theta)
S <- matrix(1e-3, n,p)

par0 <- c(Theta, M, S)

outR <- fn_optim_PLN(par0, Y, X, O)
outC <- fn_optim_Cpp(par0, Y, X, O, KY)

library(microbenchmark)
res <- microbenchmark(outR = fn_optim_PLN(par0, Y, X, O),
                      outC = fn_optim_Cpp(par0, Y, X, O, KY))
print(summary(res))
*/
