// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List fn_optim_PLNPCA_Cpp(const arma::vec par,
                               int q,
                               const arma::mat Y,
                               const arma::mat X,
                               const arma::mat O,
                               double KY) {

  int n = Y.n_rows, p = Y.n_cols, d = X.n_cols ;

  arma::mat Theta = par.subvec(0          , p*d          -1) ; Theta.reshape(p,d) ;
  arma::mat B     = par.subvec(p*d        , p*(d+q)      -1) ; B.reshape(p,q) ;
  arma::mat M     = par.subvec(p*(d+q)    , p*(d+q)+n*q  -1) ; M.reshape(n,q) ;
  arma::mat S     = par.subvec(p*(d+q)+n*q, p*(d+q)+2*n*q-1) ; S.reshape(n,q) ;

  arma::mat Z = O + X * Theta.t() + M * B.t();
//  arma::mat A = exp (Z + .5 * (S%S) * (B%B).t() ) ;
  arma::mat A = exp (Z + .5 * S * (B%B).t() ) ;

  arma::vec grd_Theta = vectorise((A-Y).t() * X);
//  arma::vec grd_B     = vectorise((A-Y).t() * M + (A.t() * (S%S)) % B) ;
  arma::vec grd_B     = vectorise((A-Y).t() * M + (A.t() * S) % B) ;
  arma::vec grd_M     = vectorise((A-Y) * B + M) ;
//  arma::vec grd_S     = vectorise(S - 1/S + (A * (B%B)) % S );
  arma::vec grd_S     = .5 * vectorise(1 - 1/S + A * (B%B) );

  arma::vec grad = join_vert(join_vert(grd_Theta, grd_B), join_vert(grd_M, grd_S)) ;

//  double objective = accu(A - Y % Z) + accu(.5 * (M%M + S%S) - log(S) - .5) + KY ;
  double objective = accu(A - Y % Z) + .5 * accu(M%M + S - log(S) - 1) + KY ;

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

## declare the objective and gradient functions for optimization
fn_optim_PLNPCA <- function(par,q,Y,X,O,KY) {

    Theta <- matrix(par[1:(p*d)]                          , p,d)
    B     <- matrix(par[p*d                + 1:(p*q)], p,q)
    M     <- matrix(par[p*(d+q)            + 1:(n*q)], n,q)
    S     <- matrix(par[p*(d+q)+(n*q) + 1:(n*q)], n,q)

    Z     <- tcrossprod(M,B) + tcrossprod(X, Theta) + O
    # A <- exp (.trunc(Z + .5*tcrossprod(S^2, B^2)))
    A <- exp (.trunc(Z + .5*tcrossprod(S, B^2)))

    gr.Theta <- crossprod(A-Y, X)
    # gr.B     <- crossprod(A-Y, M) + crossprod(A,S^2) * B
    gr.B     <- crossprod(A-Y, M) + crossprod(A,S) * B
    gr.M     <- (A-Y) %*% B + M
    # gr.S     <- S - 1/S  + (A %*% (B^2)) * S
    gr.S     <- .5 *(1 - 1/S  + .5 * A %*% (B^2))

    return(list(
#      "objective" = sum(as.numeric(A -Y*Z)) +.5*sum(M^2 + S^2 - 2*log(S)-1) + KY,
      "objective" = sum(as.numeric(A -Y*Z)) +.5*sum(M^2 + S - log(S)-1) + KY,
      "gradient"  = c(gr.Theta,gr.B,gr.M,gr.S)
    ))
  }

q <- 1
glmP  <- lapply(1:ncol(Y), function(j) glm.fit(X, Y[, j], offset = O[,j], family = poisson()))
Theta <- do.call(rbind, lapply(glmP, coefficients))
M <- matrix(runif(n*q,-1,1), n, q) ## -tcrossprod(covariates,self$init.par$Theta)
S <- matrix(runif(n*q,1e-3,1e-1), n, q)

Sigma <- cov(sapply(glmP, residuals.glm, "pearson"))
svdSigma <- svd(Sigma, nu=q, nv=0)
B <- svdSigma$u[, 1:q, drop=FALSE] %*% diag(sqrt(svdSigma$d[1:q]),nrow=q, ncol=q)

par0 <- c(Theta, B, M, S)

outR <- fn_optim_PLNPCA(par0, q, Y, X, O, KY)
outC <- fn_optim_PLNPCA_Cpp(par0, q, Y, X, O, KY)

library(microbenchmark)
res <- microbenchmark(outR = fn_optim_PLNPCA(par0, q, Y, X, O, KY),
                      outC = fn_optim_PLNPCA_Cpp(par0, q, Y, X, O, KY))
print(summary(res))
*/
