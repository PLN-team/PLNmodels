#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function \code{\link{PLNnetwork}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNnetworkfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNnetworkfamily_getModel]{getModel}} and \code{\link[=PLNnetworkfamily_plot]{plot}}. Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field penalties the sparsity level of the network in the successively fitted models
#' @field models a list of \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} object, one per penalty.
#' @field inception a \code{\link[=PLNfit-class]{PLNfit}} object, obtained when no sparsifying penalty is applied.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @field fn_optim the R functions used to compute the model's objective and gradient during the optimization process
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom nloptr nloptr
#' @importFrom glasso glasso
#' @import igraph
#' @import Matrix
#' @import ggplot2
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
PLNnetworkfamily <-
  R6Class(classname = "PLNnetworkfamily",
    inherit = PLNfamily,
     active = list(
      penalties = function() self$params
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(penalties, responses, covariates, offsets, control) {

    ## initialize fields shared by the super class
    super$initialize(responses, covariates, offsets, control)

    ## Get an appropriate grid of penalties
    if (control$trace > 0) cat("\n Recovering an appropriate grid of penalties.")
    if (is.null(penalties)) {
      range_penalties <- range(abs(self$inception$model_par$Sigma[upper.tri(self$inception$model_par$Sigma)]))
      min_log_scale <- log10(range_penalties[2])
      max_log_scale <- log10(max(range_penalties[1],range_penalties[2]*control$min.ratio))
      penalties <- 10^seq(min_log_scale, max_log_scale, len = control$nPenalties)
    } else {
      if (verbose) {
        cat("\nPenalties already set by the user")
      }
      stopifnot(all(penalties > 0))
    }

    ## instantiate as many models as penalties
    self$params <- sort(penalties, decreasing = FALSE)
    self$models <- vector("list", length(self$params))
    fit <- PLNnetworkfit$new()
    for (m in seq_along(self$params))  {
      self$models[[m]] <- fit$clone()
      self$models[[m]]$update(penalty = self$params[m])
    }

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
  Sigma <- self$inception$model_par$Sigma
  par0 <- c(self$inception$model_par$Theta,
            self$inception$var_par$M,
            self$inception$var_par$S)
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
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 0) cat("",iter)

      ## Update Omega/Sigma
      Omega <- glasso::glasso(Sigma, rho = penalty, penalize.diagonal = control$penalize.diagonal)$wi
      logDetOmega <- determinant(Omega, logarithm = TRUE)$modulus

      ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
      ## to update Theta, M and S
      opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep = "_"),
                   "maxeval"     = control$maxeval,
                   "ftol_rel"    = control$ftol_rel,
                   "ftol_abs"    = control$ftol_abs,
                   "xtol_rel"    = control$xtol_rel,
                   "xtol_abs"    = control$xtol_abs,
                   "print_level" = max(0,control$trace - 1))

      optim.out <- nloptr(par0, eval_f = private$fn_optim, lb = lower.bound, opts = opts,
                          log_detOmega = logDetOmega, Omega = Omega,
                          Y = self$responses, X = self$covariates, O = self$offsets, KY = KY)

      objective[iter]   <- optim.out$objective
      convergence[iter] <- sqrt(sum((optim.out$solution - par0)^2)/sum(par0^2))

      ## Check convergence
      if ((convergence[iter] < tol) | (iter >= maxit)) {
        cond <- TRUE
      }

      ## Post-Treatment to update Sigma
      M <- matrix(optim.out$solution[private$p*private$d             + 1:(private$n*private$p)], private$n,private$p)
      S <- matrix(optim.out$solution[(private$n+private$d)*private$p + 1:(private$n*private$p)], private$n,private$p)
      Sigma <- crossprod(M)/private$n + diag(colMeans(S), nrow = private$p, ncol = private$p)
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
    BIC <- J - (private$p * private$d + sum(Omega[upper.tri(Omega, diag = TRUE)] != 0)) * log(private$n)/2
    ICL <- BIC - .5*private$n*private$p * log(2*pi*exp(1)) - sum(log(S))

    ## Enforce symmetry of Sigma and Theta
    if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)
    if (!isSymmetric(Theta)) Omega <- Matrix::symmpart(Omega)

    self$models[[m]]$update(Omega = Omega, Sigma = Sigma, Theta = Theta,
                            M = M, S = S, J = J, BIC = BIC, ICL = ICL,
                            status = convergence[iter], iter = iter)
  }

})

PLNnetworkfamily$set("public", "optimize_approx",
  function(control) {

  ## ===========================================
  ## INITIALISATION
  ## start from the standard PLN (a.k.a. inception)
  ## keep these values for the variational parameters whatever the penalty
  ## KY <- sum(.logfactorial(self$responses)) ## constant quantity in the objective
  Theta <- self$inception$model_par$Theta
  M     <- self$inception$var_par$M
  S     <- self$inception$var_par$S
  Sigma <- crossprod(M)/private$n + diag(colMeans(S), nrow = private$p, ncol = private$p)

  for (m in seq_along(self$models))  {

    penalty <- self$models[[m]]$penalty
    ## ===========================================
    ## OPTIMISATION
    if (control$trace > 0) cat("\n sparsifying penalty =",penalty)
    if (control$trace > 1) cat("\n\t approximate version: do not optimize the variational paramters")
    if (control$trace > 1) cat("\n\t graphical-Lasso for sparse covariance estimation")

    Omega <- suppressWarnings(glasso::glasso(Sigma, rho=penalty, penalize.diagonal = control$penalize.diagonal, approx = TRUE)$wi)
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)

    ## ===========================================
    ## OUTPUT

    ## Enforce symmetry of Sigma
    if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)
    if (!isSymmetric(Theta)) Omega <- Matrix::symmpart(Omega)

    self$models[[m]]$update(Omega = Omega, Sigma = Sigma, Theta = Theta, M = M, S = S)
  }

})

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
    geom_vline(xintercept = self$getBestModel("BIC")$penalty, linetype = "dashed", alpha = 0.5)
  if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
  p
})

PLNnetworkfamily$set("public", "show",
  function() {
    super$show()
    cat(" Task: Network Inference \n")
    cat("========================================================\n")
    cat(" -", length(self$penalties) , "penalties considered: from", min(self$params), "to", max(self$params),"\n")
    cat(" - Best model (regarding BIC): lambda =", format(self$getBestModel("BIC")$penalty, digits = 3), "- R2 =", round(self$getBestModel("BIC")$criteria['R2'], 2), "\n")
  })
