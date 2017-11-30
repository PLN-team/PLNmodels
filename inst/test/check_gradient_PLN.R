rm(list=ls())
library(nloptr)
library(RcppArmadillo)
library(Rcpp)
library(PLNmodels)
source("../../R/utils.R")
source("../../R/PLN.R")

## Basic data simulation
n <- 50
p <- 20
d <- 1
rho <- 0.8
Sigma <- toeplitz(rho^(0:(p-1)))
Z <- matrix(rnorm(n*p), n,p) %*% chol(Sigma)
mu <- outer(rep(1,n), runif(p, 1, 2))
Y <- matrix(rpois(n*p, mu + exp(Z)), n,p)
X <- matrix(rep(1,n), n, 1)
O <- matrix(0,n,p)

## Set the generic
ctrl <- PLN_param(list(),n,p)
KY   <- sum(.logfactorial(Y))
X.XtXm1 <- X %*% solve(crossprod(X))
ProjOrthX <- diag(n) - tcrossprod(X.XtXm1,X)
par0 <- initializePLN(Y, X, O, ctrl) # Theta, M, S

opts <- list(
  "algorithm"   = paste("NLOPT_LD",ctrl$method, sep="_"),
  "maxeval"     = ctrl$maxeval,
  "ftol_rel"    = ctrl$ftol_rel,
  "ftol_abs"    = ctrl$ftol_abs,
  "xtol_rel"    = ctrl$xtol_rel,
  "check_derivatives" = TRUE,
  "check_derivatives_tol" = 1e-2,
  "print_level" = max(0,ctrl$trace-1)
)

cat("\n ORIGINAL PARAMETRIZATION")
x1 <- c(par0$Theta, par0$M, par0$S)
opts1 <- opts
opts1$xtol_abs <- c(rep(0, p*d), rep(0, p*n), rep(ctrl$xtol_abs, n*p))
optim1 <- nloptr(
  x1,
  eval_f = fn_optim_PLN_par1_Cpp,
  lb = c(rep(-Inf, p*d), rep(-Inf, p*n), rep(ctrl$lbvar, n*p)),
  opts = opts1,
  Y = Y, X = X, O = O, KY = KY
)

cat("\n ALTERNATIVE PARAMETRIZATION")
x2 <- c(par0$M + O + tcrossprod(X, par0$Theta), par0$S)
opts2 <- opts
opts2$xtol_abs <- c(rep(0, p*n), rep(ctrl$xtol_abs, n*p))
optim2 <- nloptr(
  x2,
  eval_f = fn_optim_PLN_par2_Cpp,
  lb     = c(rep(-Inf, p*n), rep(ctrl$lbvar, n*p)),
  opts   = opts2,
  Y = Y, ProjOrthX = ProjOrthX, O = O, KY = KY
)
cat("\n objectives:", optim1$objective, optim2$objective, "\n")

library(microbenchmark)
timings.func <- microbenchmark(
  par1 = fn_optim_PLN_par1_Cpp(x1, Y = Y, X = X, O = O, KY = KY),
  par2 = fn_optim_PLN_par2_Cpp(x2, Y = Y, ProjOrthX = ProjOrthX, O = O, KY = KY)
)
print(summary(timings.func))

opts1$check_derivatives <- FALSE
opts2$check_derivatives <- FALSE
timings.optim <- microbenchmark(
  par1 = nloptr(x1, eval_f = fn_optim_PLN_par1_Cpp, lb = c(rep(-Inf, p*d), rep(-Inf, p*n), rep(ctrl$lbvar, n*p)), opts = opts1, Y = Y, X = X, O = O, KY = KY),
  par2 = nloptr(x2, eval_f = fn_optim_PLN_par2_Cpp, lb = c(rep(-Inf, p*n), rep(ctrl$lbvar, n*p)), opts   = opts2, Y = Y, ProjOrthX = ProjOrthX, O = O, KY = KY)
)
print(summary(timings.optim))

ggplot2::autoplot(timings.func)
ggplot2::autoplot(timings.optim)

