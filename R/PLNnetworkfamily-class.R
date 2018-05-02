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
      max_pen <- max(abs(self$inception$model_par$Sigma))
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
  objective.old <- ifelse(is.na(self$inception$loglik),-Inf,-self$inception$loglik)

  ## set option for call to NLOPT: optim typ, lower bound, tolerance...
  lower.bound <- c(rep(-Inf, private$p*private$d), # Theta
                   rep(-Inf, private$n*private$p), # M
                   rep(control$lbvar, private$n*private$p)) # S
  xtol_abs <- c(rep(0, private$p*private$d), # Theta
                rep(0, private$n*private$p), # M
                rep(control$xtol_abs, private$n*private$p)) # S
  opts <- list("algorithm"   = paste("NLOPT_LD",control$method, sep = "_"),
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = xtol_abs,
               "print_level" = max(0,control$trace - 2))
  if (!control$nloptr) {
    opts$algorithm   <- control$method
    opts$lower_bound <- lower.bound
  }
  ## ===========================================
  ## GO ALONG THE PENALTY GRID (i.e the models)
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
      logDetOmega <- as.double(determinant(Omega, logarithm = TRUE)$modulus)
      if (control$nloptr) {
        optim.out <- nloptr(par0, eval_f = private$fn_optim, lb = lower.bound, opts = opts,
                            log_detOmega = logDetOmega, Omega = Omega,
                            Y = self$responses, X = self$covariates, O = self$offsets, KY = KY)
      } else {
        ## Optimize via NLOPT directly
        optim.out <- optimization_PLNnetwork(par0, self$responses, self$covariates, self$offsets, Omega, logDetOmega, opts)
        optim.out$message <- statusToMessage(optim.out$status)
      }
      objective[iter]   <- optim.out$objective + penalty * sum(abs(Omega))
      convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])

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
    ## Optimization ends with a gradient descent step rather than a glasso step.
    ## Return Sigma from glasso step to ensure that Sigma = solve(Omega)
    ## Sigma <- Sigma0 ; if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)
    ## dimnames(Sigma) <- dimnames(Omega)

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

#' Extract the regularization path of a PLNnetwork fit
#'
#' @name coefficient_path
#' @param precision a logical, should the coefficients of the precision matrix Omega or the covariance matrice Sigma be sent back. Default is \code{TRUE}.
#' @param corr a logical, should the correlation (partial in case  precision = TRUE) be sent back. Default is \code{TRUE}.
#'
#' @return  Send back a tibble/data.frame.
NULL
PLNnetworkfamily$set("public", "coefficient_path",
function(precision = TRUE, corr = TRUE) {
  lapply(self$penalties, function(x) {
    if (precision) {
      G <- self$getModel(x)$model_par$Omega
    } else {
      G <- self$getModel(x)$model_par$Sigma
      dimnames(G) <- dimnames(self$getModel(x)$model_par$Omega)
    }
    if (corr) {
      G <- ifelse(precision, -1, 1) * G / tcrossprod(sqrt(diag(G)))
    }
    G %>% melt(value.name = "Coeff", varnames = c("Node1", "Node2")) %>%
      mutate(Penalty = x,
             Node1   = as.character(Node1),
             Node2   = as.character(Node2),
             Edge    = paste0(Node1, "--", Node2)) %>%
      filter(Node1 < Node2)
  }) %>% bind_rows()
})

#' @export
PLNnetworkfamily$set("public", "density_path",
function(networks) {
  data.frame(Penalty = self$penalties,
             Density = sapply(self$penalties,
                              function(x) { self$getModel(x)$latent_network() %>% mean() }),
             stringsAsFactors = FALSE)
})

#' @export
PLNnetworkfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "EBIC"), log.x = TRUE) {
  vlines <- sapply(criteria, function(crit) self$getBestModel(crit)$penalty)
  p <- super$plot(criteria) + xlab("penalty") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
  p
})

#' @export
PLNnetworkfamily$set("public", "plot_objective",
function() {
  objective <- unlist(lapply(self$models, function(model) model$optim_par$objective))
  changes <- cumsum(unlist(lapply(self$models, function(model) model$optim_par$outer_iterations)))
  dplot <- data.frame(iteration = 1:length(objective), objective = objective)
  p <- ggplot(dplot, aes(x = iteration, y = objective)) + geom_line() +
    geom_vline(xintercept = changes, linetype="dashed", alpha = 0.25) +
    ggtitle("Objective along the alternate algorithm") + xlab("iteration (+ changes of model)") +
    annotate("text", x = changes, y = min(dplot$objective), angle = 90,
             label = paste("penalty=",format(self$criteria$param, digits = 1)), hjust = -.1, size = 3, alpha = 0.7) + theme_bw()
  p
})

PLNnetworkfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Network Inference \n")
  cat("========================================================\n")
  cat(" -", length(self$penalties) , "penalties considered: from", min(self$penalties), "to", max(self$penalties),
      "\n", "   use $penalties to see all values and access specific lambdas", "\n")
  if (!anyNA(self$criteria$EBIC))
    cat(" - Best model (regarding EBIC): lambda =", format(self$getBestModel("EBIC")$penalty, digits = 3), "- R2 =", round(self$getBestModel("EBIC")$R_squared, 2), "\n")
  if (!anyNA(self$criteria$BIC))
    cat(" - Best model (regarding BIC): lambda =", format(self$getBestModel("BIC")$penalty, digits = 3), "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  # if (!anyNA(self$criteria$ICL))
  #   cat(" - Best model (regarding ICL): lambda =", format(self$getBestModel("ICL")$penalty, digits = 3), "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNnetworkfamily$set("public", "print", function() self$show())

