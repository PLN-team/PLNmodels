#' An R6 Class to represent a collection of PLNMMfit
#'
#' @description The function \code{\link{PLNMM}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNfamily_getModel]{getModel}}, \code{\link[=plot.PLNfamily]{plot}}
#' and \code{\link[=predict.PLNfit]{predict}}.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field clusters the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNMMfit]{PLNMMfit}} object, one for each number of clusters considered.
#' @field inception a \code{\link[=PLNfit]{PLNfit}} object, obtained when one single class is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import mclust
#' @import ggplot2
#' @seealso The function \code{\link{PLNMM}}, the class \code{\link[=PLNMMfit]{PLNMMfit}}
PLNMMfamily <-
  R6Class(classname = "PLNMMfamily",
    inherit = PLNfamily,
     active = list(
      clusters = function() private$params
    )
)

PLNMMfamily$set("public", "initialize",
  function(clusters, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  private$params <- clusters

  if (control$trace > 0) cat("\n Perform GMM on the inceptive model for initializing...")

  ## instantiate as many models as clusters
  self$models <- lapply(clusters, function(k) {
    mclust_out <- mclust::Mclust(
      data           = self$inception$var_par$M,
      G              = k,
      modelNames     = "EII",
      initialization = list(hcPairs = hc(self$inception$var_par$M, "EII")),
      verbose        = FALSE)

    model <- PLNMMfit$new(
      Theta = self$inception$model_par$Theta,
      Sigma = diag(colMeans(self$inception$var_par$S)),
      Mu    = mclust_out$parameters$mean,
      M     = array(rep(self$inception$var_par$M, k), c(private$n, private$p, k)),
      S     = array(rep(self$inception$var_par$S, k), c(private$n, private$p, k)),
      Tau   = mclust_out$z
    )
    return(model)
  })

})

PLNMMfamily$set("public", "optimize",
  function(control) {

  ## set option for call to NLOPT: optim typ, lower bound, tolerance...
  opts <- list("algorithm"   = control$method,
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = c(rep(0, private$p*(private$d + 1)), # Theta + Muk
                                 rep(0, private$n*private$p), # Mk
                                 rep(control$xtol_abs, private$n*private$p)), #Sk
               "lower_bound" = c(rep(-Inf, private$p*(private$d + 1)), # Theta + Muk
                                 rep(-Inf, private$n*private$p), # Mk
                                 rep(control$lbvar, private$n*private$p)), # Sk
               "weighted"    = TRUE
  )



  ## ===========================================
  ## GO ALONG THE NUMBER OF CLUSTER (i.e the models)
  for (m in seq_along(self$models))  {
    ## k is the total number of cluster
    k <- self$models[[m]]$k

    ## ===========================================
    ## INITIALISATION
    objective_old <- -Inf
    tau_old <- matrix(NA, private$n, k); par0 <- list()
    for (k_ in 1:k) {
      tau_old[, k_] <- self$models[[m]]$posteriorProb[, k_]
      par0[[k_]] <- c(
        cbind(self$models[[m]]$model_par$Mu[,k_], self$models[[m]]$model_par$Theta),
        self$models[[m]]$var_par$M[, , k_],
        self$models[[m]]$var_par$S[, , k_]
      )
    }

    if (control$trace == 1) {
      cat("\tnumber of cluster =", k, "\r")
      flush.console()
    }
    if (control$trace > 1) {
      cat("\tnumber of cluster =", k, "- iteration:")
    }

    ## ===========================================
    ## OPTIMISATION
    cond <- FALSE; iter <- 0
    objective   <- numeric(control$maxit_out)
    objective_component <- vector("numeric", k)
    loglikObs_component <- matrix(NA, private$n, k)
    convergence <- numeric(control$maxit_out)
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 1) cat("", iter)
      ## UPDATE COMPONENTS OF THE MIXTURE
      for (k_ in 1:k) {
        ## UPDATE COMPONENTS OF THE MIXTURE
        ## weighted-PLN model - probabilities of cluster belonging
        optim_out <- optimization_PLN(
          par0[[k_]],
          self$responses,
          cbind(rep(1,private$n), self$covariates),
          self$offsets,
          tau_old[, k_],
          opts
        )
        par0[[k_]] <- optim_out$solution
        loglikObs_component[, k_] <- optim_out$loglik
      }
      ## UPDATE THE POSTERIOR PROBABILITIES
##      tau <- t(apply(sweep(loglikObs_component, 2, log(pi_hat), "+"), 1, .softmax)
      tau <- sweep(loglikObs_component, 2, log(colMeans(tau_old)), "+")
      tau <- t(apply(tau, 1, .softmax))
      tau[tau < .Machine$double.eps] <- .Machine$double.eps
      tau[tau > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
      pi_hat <- colMeans(tau)

      objective[iter] <- -sum(tau * loglikObs_component) + sum(.xlogx(tau)) - sum(tau * log(pi_hat))
      # print(- sum(tau * loglikObs_component) + sum(.xlogx(tau)) - sum(tau * log(pi)))
      # print(pi_hat)
      # print(aricode::NID(apply(tau_old, 1, which.max),apply(tau, 1, which.max)))
      convergence[iter] <- abs(objective_old - objective[iter])/abs(objective[iter])

      ## Check convergence
      if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE
      objective_old <- objective[iter]
      tau_old <- tau
    }

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output

    Theta <- list()
    Mu    <- list()
    M     <- list()
    S     <- list()
    Sigma <- list()

    for (k_ in seq.int(k)) {
      Theta[[k_]] <- matrix(par0[[k_]][1:(private$p*(private$d + 1))], private$p, private$d + 1)
      Mu[[k_]]    <- Theta[[k_]][, 1, drop = FALSE]
      Theta[[k_]] <- Theta[[k_]][, -1, drop = FALSE]
      M[[k_]]     <- matrix(par0[[k_]][private$p*(private$d + 1)   + 1:(private$n*private$p)], private$n,private$p)
      S[[k_]]     <- matrix(par0[[k_]][private$p*(private$d + 1 + private$n) + 1:(private$n*private$p)], private$n, private$p)
      Sigma[[k_]] <- t(M[[k_]]) %*% diag(tau[, k_]) %*% M[[k_]]/ sum(tau[, k_]) + diag(tau[, k_] %*% S[[k_]])
    }

    Theta <- Reduce("+", Theta)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    Sigma <- Reduce("+", Sigma)
    rownames(Sigma) <- colnames(Sigma) <- colnames(self$responses)
    M <- array(unlist(M), c(private$n, private$p, k) )
    S <- array(unlist(S), c(private$n, private$p, k) )

    self$models[[m]]$update(
      Sigma = Sigma, Theta = Theta, M = M, S = S, J = -objective[iter], Tau = tau,
          monitoring = list(objective        = objective[1:iter],
                            convergence      = convergence[1:iter],
                            outer_iterations = iter))
    if (control$trace > 1) {
      cat("\r                                                                                    \r")
      flush.console()
    }

  }
})

PLNMMfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "ICL")) {
  vlines <- sapply(criteria, function(crit) self$getBestModel(crit)$k)
  p <- super$plot(criteria) + xlab("# of clusters") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  p
})

PLNMMfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Mixture Model \n")
  cat("========================================================\n")
  cat(" - Number of clusters considered: from", min(self$clusters), "to", max(self$clusters),"\n")
  cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$cluster, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$cluster, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNMMfamily$set("public", "print", function() self$show())
