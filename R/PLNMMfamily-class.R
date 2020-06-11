#' An R6 Class to represent a collection of PLNMMfit
#'
#' @description The function [PLNMM()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()], [getModel()] and [`plot()`][plot.PLNPCAfamily()].
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param control a list for controlling the optimization. See details.
#' @param clusters the dimensions of the successively fitted models
#' @param model model used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call #'
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

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation -----------------------
    #' @description Initialize all models in the collection.
      initialize = function(clusters, responses, covariates, offsets, model, xlevels, control) {

        ## initialize the required fields
        super$initialize(responses, covariates, offsets, rep(1, nrow(responses)), control)
        private$params <- clusters

        if (control$trace > 0) cat("\n Perform GMM on the inceptive model for initializing...")
        myPLN <- PLNfit$new(responses, covariates, offsets, rep(1, nrow(responses)), model, xlevels, control)
        myPLN$optimize(responses, covariates, offsets, rep(1, nrow(responses)), control)

        data0 <- myPLN$var_par$M + tcrossprod(covariates, myPLN$model_par$Theta)
        initMclust <- mclust::hc(data0, "EII")

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
          PLNMMfit$new(inception = myPLN, tau = mclust_out$z)
        })
      },
      optimize = function(control) {
        ## ===========================================
        ## GO ALONG THE NUMBER OF CLUSTER (i.e the models)
        for (m in seq_along(self$models))  {
          ## k is the total number of cluster
          k <- self$models[[m]]$k

          ## ===========================================
          ## INITIALISATION
          tau <- self$models[[m]]$posteriorProb
          loglikObs_component <- sapply(self$models[[m]]$components, function(model) model$loglik_vec)
          objective    <- numeric(control$maxit_out)
          objective[1] <- -sum(tau * loglikObs_component) + sum(.xlogx(tau )) - sum(tau * log(self$models[[m]]$mixtureParam))
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

            ## UPDATE THE MIXTURE MODEL VIA OPTIMIZATION OF PLNMM
            self$models[[m]]$optimize(self$responses, self$covariates, self$offsets, tau, control)
            loglikObs_component <- sapply(self$models[[m]]$components, function(model) model$loglik_vec)

            ## UPDATE THE POSTERIOR PROBABILITIES
            if (k > 1) { # only needed when at least 2 components!
              tau <- t(apply(sweep(loglikObs_component, 2, colMeans(tau), "+"), 1, .softmax))
              tau[tau < .Machine$double.eps] <- .Machine$double.eps
              tau[tau > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
            }
            objective[iter + 1] <- -sum(tau * loglikObs_component) + sum(.xlogx(tau)) - sum(tau * log(colMeans(tau)))
            convergence[iter]   <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter])

            ## Check convergence
            if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE
          }

          ## ===========================================
          ## OUTPUT
          ## formating parameters for output

          self$models[[m]]$update(tau = tau, J = -objective[iter + 1],
            monitoring = list(objective        = objective[1:(iter + 1)],
                              convergence      = convergence[1:iter],
                              outer_iterations = iter))
          if (control$trace > 1) {
            cat("\r                                                                                    \r")
            flush.console()
          }

        }
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -------------
      #' @description
      #' Lineplot of selected criteria for all models in the collection
      #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
      #' @param annotate Logical. Should R2 be added to the plot (defaults to `TRUE`)
      #' @return A [`ggplot2`] object
      plot = function(criteria = c("loglik", "BIC", "ICL"), annotate = FALSE) {
        # vlines <- sapply(intersect(criteria, c("BIC", "ICL")), function(crit) self$getBestModel(crit)$k)
        p <- super$plot(criteria, annotate) + xlab("# of clusters") # + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
        p
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Extractors   -------------------
      #' @description Extract best model in the collection
      #' @param crit a character for the criterion used to performed the selection. Either
      #' "BIC", "ICL", or "R_squared". Default is `BIC`
      #' @return a [`PLNPCAfit`] object
      getBestModel = function(crit = c("BIC", "ICL", "R_squared")){
        crit <- match.arg(crit)
        stopifnot(!anyNA(self$criteria[[crit]]))
        id <- 1
        if (length(self$criteria[[crit]]) > 1) {
          id <- which.max(self$criteria[[crit]])
        }
        model <- self$models[[id]]$clone()
        model
      },
      show = function() {
        super$show()
        cat(" Task: Mixture Model \n")
        cat("========================================================\n")
        cat(" - Number of clusters considered: from", min(self$clusters), "to", max(self$clusters),"\n")
        cat(" - Best model (regarding BIC): cluster =", self$getBestModel("BIC")$k, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
        cat(" - Best model (regarding ICL): cluster =", self$getBestModel("ICL")$k, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
      },
      print = function() self$show()
    )
)

