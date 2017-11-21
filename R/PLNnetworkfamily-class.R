#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function \code{\link{PLNnetwork}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNnetworkfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNnetworkfamily_getModel]{getModel}} and \code{\link[=PLNnetworkfamily_plot]{plot}}.
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
      penalties = function() private$params
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(penalties, responses, covariates, offsets, control) {

    ## initialize fields shared by the super class
    super$initialize(responses, covariates, offsets, control)

    ## Get an appropriate grid of penalties
    if (is.null(penalties)) {
      if (control$trace > 1) cat("\n Recovering an appropriate grid of penalties.")
      max_pen <- max(abs(self$inception$model_par$Sigma[upper.tri(self$inception$model_par$Sigma)]))
      penalties <- 10^seq(log10(max_pen), log10(max_pen*control$min.ratio), len = control$nPenalties)
    } else {
      if (control$trace > 1) cat("\nPenalties already set by the user")
      stopifnot(all(penalties > 0))
    }

    ## instantiate as many models as penalties
    private$params <- sort(penalties, decreasing = TRUE)
    self$models <- lapply(private$params, function(penalty) {
      PLNnetworkfit$new(penalty = penalty)
    })

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

  ## start from the standard PLN (a.k.a. inception)
  Sigma <- self$inception$model_par$Sigma
  par0  <- c(self$inception$model_par$Theta,
             self$inception$var_par$M,
             self$inception$var_par$S)
  Omega  <- diag(1/diag(Sigma), nrow = private$p, ncol = private$p)
  Sigma0 <- diag(  diag(Sigma), nrow = private$p, ncol = private$p)
  objective.old <- -self$inception$loglik

  ## set option for call to NLOPT: optim typ, lower bound, tolerance...
  lower.bound <- c(rep(-Inf, private$p*private$d), # Theta
                   rep(-Inf, private$n*private$p), # M
                   rep(control$lbvar, private$n*private$p)) # S
  xtol_abs <- c(rep(0, private$p*private$d), # Theta
                rep(0, private$n*private$p), # M
                rep(control$xtol, private$n*private$p)) # S
  opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep = "_"),
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = xtol_abs,
               # "check_derivatives" = TRUE,
               "print_level" = max(0,control$trace - 2))

  ## ===========================================
  ## GET ALONG THE PENLATY GRID (i.e the models)

  for (m in seq_along(self$models))  {

    penalty <- self$models[[m]]$penalty

    ## ===========================================
    ## OPTIMISATION
    if (control$trace == 1) {
      cat("\tsparsifying penalty =",penalty, "\r")
      flush.console()
    }
    if (control$trace > 1) {
      cat("\tsparsifying penalty =",penalty, "- iteration:")
    }

    cond <- FALSE; iter <- 0
    objective   <- numeric(control$maxit_out)
    convergence <- numeric(control$maxit_out)
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 1) cat("",iter)

      ## CALL TO GLASSO TO UPDATE Omega/Sigma
      glasso_out <- suppressWarnings(
                      glasso(Sigma,
                             rho = penalty,
                             penalize.diagonal = control$penalize.diagonal,
                             start = ifelse(control$warm, "warm", "cold"), w.init = Sigma0, wi.init = Omega
                             )
                      )
      Omega  <- glasso_out$wi ; if (!isSymmetric(Omega)) Omega <- Matrix::symmpart(Omega)
      Sigma0 <- glasso_out$w

      ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
      logDetOmega <- determinant(Omega, logarithm = TRUE)$modulus
      optim.out <- nloptr(par0, eval_f = private$fn_optim, lb = lower.bound, opts = opts,
                          log_detOmega = logDetOmega, Omega = Omega,
                          Y = self$responses, X = self$covariates, O = self$offsets, KY = KY)
      objective[iter]   <- optim.out$objective + penalty * sum(abs(Omega))
      convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective.old)

      ## Check convergence
      if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE

      ## Post-Treatment to update Sigma
      M <- matrix(optim.out$solution[private$p*private$d               + 1:(private$n*private$p)], private$n,private$p)
      S <- matrix(optim.out$solution[(private$n + private$d)*private$p + 1:(private$n*private$p)], private$n,private$p)
      Sigma <- crossprod(M)/private$n + diag(colMeans(S), nrow = private$p, ncol = private$p)
      par0 <- optim.out$solution
      objective.old <- objective[iter]
    }

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(private$p*private$d)], private$p, private$d)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    dimnames(S)     <- dimnames(self$responses)
    dimnames(M)     <- dimnames(self$responses)
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)

    # ## Enforce symmetry of Sigma and Theta
    # if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)

    self$models[[m]]$update(Omega = Omega, Sigma = Sigma, Theta = Theta, M = M, S = S, J = -optim.out$objective,
                            monitoring = list(objective = objective[1:iter],
                                              convergence = convergence[1:iter],
                                              outer_iterations = iter,
                                              inner_iterations = optim.out$iterations,
                                              inner_status = optim.out$status,
                                              inner_message = optim.out$message))
    if (control$trace > 1) {
      cat("\r                                                                                    \r")
      flush.console()
    }

  }

})

PLNnetworkfamily$set("public", "optimize_MB",
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
    if (control$trace > 0) {
      cat("\tsparsifying penalty =",penalty, "\r")
      flush.console()
    }

    Omega <- suppressWarnings(glasso::glasso(Sigma, rho = penalty, penalize.diagonal = control$penalize.diagonal, approx = TRUE)$wi)
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)

    ## ===========================================
    ## OUTPUT

    ## Enforce symmetry of Sigma
    if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)
    if (!isSymmetric(Theta)) Omega <- Matrix::symmpart(Omega)

    self$models[[m]]$update(Omega = Omega, Sigma = Sigma, Theta = Theta, M = M, S = S,
                            monitoring = list(objective = NA))
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
  stopifnot(!anyNA(self$criteria))
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
  cat(" -", length(self$penalties) , "penalties considered: from", min(self$penalties), "to", max(self$penalties),"\n")
  if (!anyNA(self$criteria$BIC))
    cat(" - Best model (regarding BIC): lambda =", format(self$getBestModel("BIC")$penalty, digits = 3), "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  if (!anyNA(self$criteria$ICL))
    cat(" - Best model (regarding ICL): lambda =", format(self$getBestModel("ICL")$penalty, digits = 3), "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNnetworkfamily$set("public", "print", function() self$show())
