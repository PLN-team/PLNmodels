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
#' @param model model used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call and used for predictions.
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#'
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
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
    initialize = function(ranks, responses, covariates, offsets, weights, model, xlevels, control) {
      ## initialize the required fields
      super$initialize(responses, covariates, offsets, weights, control)
      private$params <- ranks

      ## save some time by using a common SVD to define the inceptive models
      M <- do.call(cbind, lapply(1:ncol(responses), function(j)
        residuals(lm.wfit(covariates, log(1 + responses[,j]), w = weights, offset = offsets[, j]))))
      control$svdM <- svd(M, nu = max(ranks), nv = ncol(responses))

      ## instantiate as many models as ranks
      self$models <- lapply(ranks, function(rank){
        model <- PLNPCAfit$new(rank, responses, covariates, offsets, weights, model, xlevels, control)
        model
      })
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization -------------------
    #' @description Call to the C++ optimizer on all models of the collection
    optimize = function(control) {
      self$models <- mclapply(self$models, function(model) {
        if (control$trace == 1) {
          cat("\t Rank approximation =",model$rank, "\r")
          flush.console()
        }
        if (control$trace > 1) {
          cat(" Rank approximation =",model$rank)
          cat("\n\t conservative convex separable approximation for gradient descent")
        }
        model$optimize(self$responses, self$covariates, self$offsets, self$weights, control)
        model
      }, mc.cores = control$cores, mc.allow.recursive = FALSE)
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

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methdods -------------
    #' @description
    #' Lineplot of selected criteria for all models in the collection
    #' @param criteria A valid model selection criteria for the collection of models. Any of "loglik", "BIC" or "ICL" (all).
    #' @param annotate Logical. Should R2 be added to the plot (defaults to `TRUE`)
    #' @return A [`ggplot2`] object
    plot = function(criteria = c("loglik", "BIC", "ICL"), annotate = TRUE) {
      vlines <- sapply(intersect(criteria, c("BIC", "ICL")) , function(crit) self$getBestModel(crit)$rank)
      p <- super$plot(criteria, annotate) + xlab("rank") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
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
      cat(" - Best model (greater BIC): rank = ", self$getBestModel("BIC")$rank, " - R2 = ", round(self$getBestModel("BIC")$R_squared, 2), "\n", sep = "")
      cat(" - Best model (greater ICL): rank = ", self$getBestModel("ICL")$rank, " - R2 = ", round(self$getBestModel("ICL")$R_squared, 2), "\n", sep = "")
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
