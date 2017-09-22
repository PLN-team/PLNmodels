#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function \code{\link{PLNnetwork}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNnetworkfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNnetworkfamily_getModel]{getModel}} and \code{\link[=PLNnetworkfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' @field penalties the sparsity level of the network in the successively fitted models
#' @field models a list of \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} object, one per penalty.
#' @field inception a \code{\link[=PLNfit-class]{PLNfit}} object, obtained when no sparsifying penalty is applied.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field fn_optim the R functions used to compute the model's objective and gradient during the optimization process
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom nloptr nloptr
#' @importFrom glasso glasso
#' @import ggplot2
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
PLNnetworkfamily <-
  R6Class(classname = "PLNnetworkfamily",
    inherit = PLNfamily,
     public = list(
      penalties = "numeric"
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(nModels, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  KY <-sum(.logfactorial(self$responses))

  ## instantiate as many models as penalties
  fit <- PLNnetworkfit$new()
  self$models <- vector("list", nModels)
  for (m in 1:nModels) {
    self$models[[m]] <- fit$clone()
  }
  # names(self$models) <- as.character(round())

  ## declare the objective and gradient functions for optimization
  private$fn_optim <- fn_optim_PLNnetwork_Cpp
})

# One (coordinate-descent) optimization for each \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} model, using
# the same starting point (inception model from \code{\link[=PLNfit-class]{PLNfit}}) for all models.
#
PLNnetworkfamily$set("public", "optimize",
  function(control) {

  ## ===========================================
  ## INITIALISATION
  KY <- sum(.logfactorial(self$responses)) ## constant quantity in the objective
  maxit <- control$out.maxit
  tol <- control$out.tol

  ## INITIALIZATION: start from the standard PLN (a.k.a. inception)
  Sigma <- self$inception$model.par$Sigma
  par0 <- c(self$inception$model.par$Theta,
            self$inception$variational.par$M,
            self$inception$variational.par$S)
  lower.bound <- c(rep(-Inf, private$p*private$d), # Theta
                   rep(-Inf, private$n*private$p), # M
                   rep(control$lbvar, private$n*private$p)) # S

  for (m in seq_along(self$models))  {

    penalty <- self$models[[m]]$penalty

    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n sparsifying penalty =",penalty)
    if (control$trace > 1) cat("\n\t conservative convex separable approximation for gradient descent")
    if (control$trace > 1) cat("\n\t graphical-Lasso for sparse covariance estimation")

    cond <- FALSE; iter <- 0
    convergence <- numeric(maxit)
    objective   <- numeric(maxit)
    if (control$trace > 0) cat("\n\titeration: ")
    while(!cond) {
      iter <- iter + 1
      if (control$trace > 0) cat("",iter)

      ## Update Omega/Sigma
      Omega <- glasso::glasso(Sigma, rho=penalty, penalize.diagonal = control$penalize.diagonal)$wi
      logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus

      ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
      ## to update Theta, M and S
      opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep="_"),
                   "maxeval"     = control$maxeval,
                   "ftol_rel"    = control$ftol_rel,
                   "ftol_abs"    = control$ftol_abs,
                   "xtol_rel"    = control$xtol_rel,
                   "xtol_abs"    = control$xtol_abs,
                   "print_level" = max(0,control$trace-1))

      optim.out <- nloptr(par0, eval_f = private$fn_optim, lb = lower.bound, opts = opts,
                          log_detOmega=logDetOmega, Omega=Omega,
                          Y=self$responses, X=self$covariates, O=self$offsets, KY=KY)

      objective[iter]   <- optim.out$objective
      convergence[iter] <- sqrt(sum((optim.out$solution - par0)^2)/sum(par0^2))

      ## Check convergence
      if ((convergence[iter] < tol) | (iter > maxit)) {
        cond <- TRUE
      }

      ## Post-Treatment to update Sigma
      M <- matrix(optim.out$solution[private$p*private$d          + 1:(private$n*private$p)], private$n,private$p)
      S <- matrix(optim.out$solution[(private$n+private$d)*private$p + 1:(private$n*private$p)], private$n,private$p)
      Sigma <- crossprod(M)/private$n + diag(colMeans(S))
      par0 <- optim.out$solution

    }
    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(private$p*private$d)], private$p, private$d)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    dimnames(S)     <- dimnames(self$responses)
    dimnames(M)     <- dimnames(self$responses)
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)

    ## compute some criteria for evaluation
    J   <- -optim.out$objective
    BIC <- J - .5 * (private$p * private$d + sum(Omega[upper.tri(Omega, diag = FALSE)]!=0)) * log(private$n)
    ICL <- BIC - .5*private$n*private$p *log(2*pi*exp(1)) - sum(log(S))

    self$models[[m]]$model.par       <- list(Omega = Omega, Sigma = Sigma, Theta = Theta)
    self$models[[m]]$variational.par <- list(M = M, S = S)
    self$models[[m]]$criteria        <- c(J = J, BIC = BIC, ICL = ICL)
    self$models[[m]]$convergence     <- data.frame(convergence = convergence[1:iter], objective = objective[1:iter])
  }

})

#' One fit for each \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} model, using
#' the same starting point (inception model from \code{\link[=PLNfit-class]{PLNfit}}) for all models.
#' Unlike \code{\link[=PLNnetworkfamily_optimize]{optimize}}, the optimization does not use an iterative procedure:
#' - Theta, M and S are fixed to their inception values
#' - Omega/Sigma is optimized only once using graphical lasso
#'
#' @name PLNnetworkfamily_optimize_approx
#'
#' @importFrom glasso glasso
PLNnetworkfamily$set("public", "optimize_approx",
  function(control) {

  ## ===========================================
  ## INITIALISATION
  ## start from the standard PLN (a.k.a. inception)
  ## keep these values for the variational parameters whatever the penalty
  KY <- sum(.logfactorial(self$responses)) ## constant quantity in the objective
  Theta <- self$inception$model.par$Theta
  M     <- self$inception$variational.par$M
  S     <- self$inception$variational.par$S
  Sigma <- crossprod(M)/private$n + diag(colMeans(S))
  par0 <- c(Theta, M, S)

  for (m in seq_along(self$models))  {

    penalty <- self$models[[m]]$penalty
    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n sparsifying penalty =",penalty)
    if (control$trace > 1) cat("\n\t approximate version: do not optimize the variational paramters")
    if (control$trace > 1) cat("\n\t graphical-Lasso for sparse covariance estimation")
    Omega <- glasso::glasso(Sigma, rho=penalty, penalize.diagonal = control$penalize.diagonal)$wi
    logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)

    ## ===========================================
    ## OUTPUT
    ## compute some criteria for evaluation
    J   <- -self$fn_optim(par0, logDetOmega, Omega, self$responses, self$covariates, self$offsets, KY)$objective
    BIC <- J - (private$p * private$d + .5*sum(Omega[upper.tri(Omega, diag = FALSE)]!=0)) * log(private$n)
    ICL <- BIC - .5*private$n*private$p *log(2*pi*exp(1)) - .5*sum(log(S))

    self$models[[m]]$model.par       <- list(Omega = Omega, Sigma = Sigma, Theta = Theta)
    self$models[[m]]$variational.par <- list(M = M, S = S)
    self$models[[m]]$criteria        <- c(J = J, BIC = BIC, ICL = ICL)
    self$models[[m]]$convergence     <- data.frame(convergence = NA, objective = -J)
  }

})

## JC: I don't want a documented function in a separated file for that

# Set penalties in a \code{\link[=PLNnetworkfamily-class]{PLNnetworkfamily}} family
# of \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} models
#
#
# param penalties A numeric vector with the sparsity levels of the networks in family. If NULL, a relevant
#                  vector is automatically computed from the inception model.
# param nPenalties The number of penalties to use. A warning is thrown if this does not match the number of
#                   models in the family.
# param verbose Logical. Controls the amount of screen output.

PLNnetworkfamily$set("public", "setPenalties",
  function(penalties, nPenalties, verbose) {
    if (is.null(penalties)) {
      cov.unpenalized <- self$inception$model.par$Sigma
      range.penalties <- range(abs(cov.unpenalized[upper.tri(cov.unpenalized)]))
      penalties <- 10^seq(log10(range.penalties[1]), log10(range.penalties[2]), len=nPenalties)
    } else {
      if (verbose) {
        cat("\nPenalties already set by the user")
      }
      nPenalties <- length(penalties)
      stopifnot(all(penalties > 0))
    }
    ## Compare number of penalties and number of models in family
    if (nPenalties != length(self$models)) {
      warning(paste0("The number of penalties (", nPenalties, ") does not match the number of models (",
                     length(self$models), ") in the family."))
    }
    ## sort penalties and round
    self$penalties <- round(sort(penalties, decreasing = FALSE),16)
    names(self$models) <- as.character(self$penalties)
    for (m in seq_along(self$models))  {
      self$models[[m]]$penalty <- self$penalties[m]
    }
    invisible(self)
  }
)

#' A plot method for a collection of PLNnetworkfit
#'
#' @name PLNnetworkfamily_plot
#' @import ggplot2
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC.
NULL
PLNnetworkfamily$set("public", "plot",
function(log.x=TRUE) {
  p <- super$plot() + xlab("penalty") +
    geom_vline(xintercept=self$getBestModel("BIC")$penalty, linetype="dashed", alpha=0.5)
  if (log.x) p <- p + ggplot2::coord_trans(x="log10")
  p
})

PLNnetworkfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Network Inference\n")
  cat("======================================================\n")
  cat(paste(" -", length(public$penalties) , "penalties considered: from", format(min(public$penalties),digits=3), "to", format(max(public$penalties),digits=3),"\n"))
  cat(" - Best model (regardings BIC): penalty =", format(self$getBestModel("BIC")$penalty,digits=3), "\n")
})
