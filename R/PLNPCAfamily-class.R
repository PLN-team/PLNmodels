#' An R6 Class to represent a collection of PLNPCAfit
#'
#' @description The function [PLNPCA()] produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for [getBestModel()],
#' [getModel()] and [`plot()`][plot.PLNPCAfamily()].
#'
## Parameters shared by many methods
#' @param ranks the dimensions of the successively fitted models
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights the vector of observation weights
#' @param control list controlling the optimization and the model
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
#' @import future
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#' class(myPCAs)
#' @seealso The function [PLNPCA()], the class [PLNPCAfit()]
PLNPCAfamily <- R6Class(
  classname = "PLNPCAfamily",
  inherit = PLNfamily,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation -----------------------
    #' @description Initialize all models in the collection.
    initialize = function(ranks, responses, covariates, offsets, weights, formula, control) {
      ## initialize the required fields
      super$initialize(responses, covariates, offsets, weights, control)
      private$params <- ranks
      ## save some time by using a common SVD to define the inceptive models
      control$inception <- PLNfit$new(responses, covariates, offsets, weights, formula, control)
      control$svdM <- svd(control$inception$var_par$M, nu = max(ranks), nv = ncol(responses))
      ## instantiate as many models as ranks
      self$models <- lapply(ranks, function(rank){
        model <- PLNPCAfit$new(rank, responses, covariates, offsets, weights, formula, control)
        model
      })
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization -------------------
    #' @description Call to the C++ optimizer on all models of the collection
    #' @param config list controlling the optimization.
    optimize = function(config) {
      self$models <- future.apply::future_lapply(self$models, function(model) {
        if (config$trace == 1) {
          cat("\t Rank approximation =",model$rank, "\r")
          flush.console()
        }
        if (config$trace > 1) {
          cat(" Rank approximation =",model$rank)
          cat("\n\t conservative convex separable approximation for gradient descent")
        }
        model$optimize(self$responses, self$covariates, self$offsets, self$weights, config)
        model
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors   -------------------
    #' @description Extract model from collection and add "PCA" class for compatibility with [`factoextra::fviz()`]
    # @inheritParams getModel
    #' @param var	value of the parameter (rank for PLNPCA, sparsity for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
    #' @param index Integer index of the model to be returned. Only the first value is taken into account.
    #' @return a [`PLNPCAfit`] object
    getModel = function(var, index = NULL) {
      model <- super$getModel(var, index)
      class(model) <- c(class(model)[class(model) != "R6"], "PCA", "R6")
      model
    },
    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection. Either
    #' "ICL", "BIC". Default is `ICL`
    #' @return a [`PLNPCAfit`] object
    getBestModel = function(crit = c("ICL", "BIC")){
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id <- which.max(self$criteria[[crit]])
      }
      self$getModel(index = id)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods -------------
    #' @description
    #' Lineplot of selected criteria for all models in the collection
    #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
    #' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
    #' (loglik - penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the natural direction, on
    #' the same scale as the log-likelihood.
    #' @return A [`ggplot2::ggplot`] object
    plot = function(criteria = c("loglik", "BIC", "ICL"), reverse = FALSE) {
      vlines <- sapply(intersect(criteria, c("BIC", "ICL")) , function(crit) self$getBestModel(crit)$rank)
      p <- super$plot(criteria, reverse) + xlab("rank") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      p
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print methods ------------------
    #' @description User friendly print method
    show = function() {
      super$show()
      cat(" Task: Principal Component Analysis\n")
      cat("========================================================\n")
      cat(" - Ranks considered: from ", min(self$ranks), " to ", max(self$ranks),"\n", sep = "")
      cat(" - Best model (greater BIC): rank = ", self$getBestModel("BIC")$rank, "\n", sep = "")
      cat(" - Best model (greater ICL): rank = ", self$getBestModel("ICL")$rank, "\n", sep = "")
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## End of methods -----------------

  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field ranks the dimensions of the successively fitted models
    ranks = function() private$params
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
