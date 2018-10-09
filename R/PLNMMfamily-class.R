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
     ),
    public  = list(
      initialize = function(clusters, responses, covariates, offsets, control) {

        control$inception <- "PLN"
        super$initialize(responses, covariates, offsets, control)
        private$params <- clusters

        if (control$trace > 0) cat("\n Perform GMM on the inceptive model for initializing...")

        data0 <- self$inception$var_par$M + tcrossprod(covariates, self$inception$model_par$Theta)
        initMclust <- hc(data0, "EII")

        ## instantiate as many PLNMMfit as choices for the number of components
        self$models <- lapply(clusters, function(k) {
          mclust_out <- mclust::Mclust(
            data           = data0,
            G              = k,
            modelNames     = "EII",
            initialization = list(hcPairs = initMclust),
            verbose        = FALSE)
          ## each PLNMMfit will itself instantiate as many PLNmodels
          ## as the current choice of number of components
          PLNMMfit$new(inception = self$inception, tau = mclust_out$z)
        })
      }
    )
)

PLNMMfamily$set("public", "optimize",
  function(control) {

  ## set option for call to NLOPT: optim type, lower bound, tolerance...
  opts <- list("algorithm"   = control$method,
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = c(rep(0, private$p*(private$d)), # Theta_k
                                 rep(0, private$n*private$p),   # Mk
                                 rep(control$xtol_abs, private$n*private$p)), #Sk
               "lower_bound" = c(rep(-Inf, private$p*(private$d)), # Theta_k
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
    tau_old <- self$models[[m]]$posteriorProb
    loglikObs_component <- sapply(self$models[[m]]$components, function(model) model$loglik_vec)
    objective    <- numeric(control$maxit_out)
    objective[1] <- -sum(tau_old * loglikObs_component) + sum(.xlogx(tau_old )) - sum(tau_old * log(self$models[[m]]$mixtureParam))
    convergence <- numeric(control$maxit_out)

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
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 1) cat("", iter)

      ## UPDATE COMPONENTS OF THE MIXTURE
      ## weighted-PLN models- with weights equal to probabilities of cluster belonging
      self$models[[m]]$optimize(self$responses, self$covariates, self$offsets, tau_old, opts)
      loglikObs_component <- sapply(self$models[[m]]$components, function(model) model$loglik_vec)
      ## Stop here for a single class

      ## UPDATE THE POSTERIOR PROBABILITIES
      if (k > 1) {
        tau <- t(apply(sweep(loglikObs_component, 2, colMeans(tau_old), "+"), 1, .softmax))
        tau[tau < .Machine$double.eps] <- .Machine$double.eps
        tau[tau > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
      } else {
        tau <- tau_old
      }
      objective[iter + 1] <- -sum(tau * loglikObs_component) + sum(.xlogx(tau)) - sum(tau * log(colMeans(tau)))

      # print(- sum(tau * loglikObs_component) + sum(.xlogx(tau)) - sum(tau * log(pi)))
      # print(pi_hat)
      # print(aricode::NID(apply(tau_old, 1, which.max),apply(tau, 1, which.max)))
      convergence[iter] <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter])

      ## Check convergence
      if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE
      objective_old <- objective[iter]
      tau_old <- tau
    }

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output

    self$models[[m]]$update(tau = tau, J = -objective[iter+1],
      monitoring = list(objective        = objective[1:iter+1],
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
