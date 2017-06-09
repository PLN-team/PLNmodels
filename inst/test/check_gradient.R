rm(list=ls())
load("~/svn/sparsepca/Pgm/PCA/Output/Corinne_data.RData")
source("~/git/PLNmodels/R/utils.R")
library(nloptr)

formula <- Data$count ~ 1 + offset(log(Data$offset))

frame  <- model.frame(formula)
Y      <- model.response(frame)
X      <- model.matrix(formula)
O      <- model.offset(frame)
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

KY <- sum(.logfactorial(Y)) ## constant quantity in the objective

func_optim <- function(par) {
  Theta <- matrix(par[1:(p*d)]                         , p,d)
  M     <- matrix(par[p*d          + 1:(n*p)], n,p)
  S     <- matrix(par[(n+d)*p + 1:(n*p)], n,p)

  Omega <- solve(crossprod(M)/n + diag(colMeans(S)))
  logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus

  Z <- O + tcrossprod(X, Theta) + M
  A <- exp (.trunc(Z + .5*S))

  logP.Z  <- n/2 * (logDetOmega - sum(diag(Omega)*colMeans(S))) - .5*sum(diag(Omega %*% crossprod(M)))

  obj <- sum(as.numeric(A - Y*Z)) - logP.Z - .5*sum(log(S)+1) + KY

  gr.Theta <- crossprod(X, A - Y)
  gr.M  <- M %*% Omega + A - Y
  gr.S  <- .5 * (matrix(rep(diag(Omega),n), n, p, byrow = TRUE) + A - 1/S)

  grad <- c(gr.Theta, gr.M, gr.S)

  return(list("objective" = obj, "gradient" = grad))
}

glmP  <- lapply(1:ncol(Y), function(j) glm.fit(X, Y[, j], offset = O[,j], family = poisson()))
Theta <- do.call(rbind, lapply(glmP, coefficients))
M <- matrix(0, n, p) ## -tcrossprod(covariates,self$init.par$Theta)
S <- matrix(1e-3, n,p)

par0 <- c(Theta, M, S)

lower.bound <- c(rep(-Inf, p*d), # Theta
                 rep(-Inf, n*p), # M
                 rep(1e-5, n*p)) # S

upper.bound <- rep(Inf, p*d+2*n*p)

# opts$local_opts <- list("algorithm" = "NLOPT_LD_LBFGS",
#                         xtol_rel = 1e-6,  eval_grad_f = gradient.unpenalized)


# check <- check.derivatives(par0,
#                   objective.unpenalized,
#                   gradient.unpenalized,
#                   check_derivatives_tol = 1e-2,
#                   KY = KY)

## LBFGS VAR1, VAR2 ok MMA suoer rapide mais fragile
opts <- list("algorithm" = NULL,
             "maxeval"   = 10000,
             "xtol_rel"  = 1e-4,
             "ftol_abs"  = 1e-6,
             "ftol_rel"  = 1e-6,
             "print_level" = 0)

# opts$algorithm = "NLOPT_LD_LBFGS"
# lbfgs.out <- nloptr(par0,
#                       eval_f      = objective.unpenalized,
#                       eval_grad_f = gradient.unpenalized,
#                       lb = lower.bound, ub = upper.bound, opts = opts)

opts$algorithm = "NLOPT_LD_MMA"
t <- system.time(mma.out <- nloptr(par0,
                 eval_f      = func_optim,
                 lb = lower.bound, ub = upper.bound, opts = opts))

# opts$algorithm = "NLOPT_LD_VAR1"
# var1.out <- nloptr(par0,
#                       eval_f      = objective.unpenalized,
#                       eval_grad_f = gradient.unpenalized,
#                       lb = lower.bound, ub = upper.bound, opts = opts)
#
# opts$algorithm = "NLOPT_LD_VAR2"
# var2.out <- nloptr(par0,
#                       eval_f      = objective.unpenalized,
#                       eval_grad_f = gradient.unpenalized,
#                       lb = lower.bound, ub = upper.bound, opts = opts)
#
