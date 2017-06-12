#' An R6 Class to represent a PLNfit in a network inference framework
#'
#' @description The function \code{\link{PLNnetwork}} produces a collection of models which are instances of object with class \code{PLNnetworkfit}.
#'
#' This class (will) come with a set of methods, some of them being useful for the user
#'
#' Fields should not be changed or manipulated by the user as they are updated internally.
#'
#' @field penalty the level of sparsity in the current model
#' @field model.par a list with two matrices, B and Theta, which are the estimated parameters of the pPCA model
#' @field variation.par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field criteria a named vector with the value of some criteria (variational lower bound J, BIC, ICL, R2, lmin and lmax) for the different models.
#' @field convergence numeric; the convergence status of the L-BFGS-B method
#' @field loglik numeric; the Poisson loglikelihood of the current model
#' @include PLNnetworkfit-class.R
#' @importFrom R6 R6Class
#' @importFrom glasso glasso
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}}
PLNnetworkfit <-
  R6Class(classname = "PLNnetworkfit",
    inherit = PLNfit,
    public  = list(
      penalty    = NULL,
      initialize = function(penalty = NA,
                            model.par=NA, variational.par=NA, criteria=NA, convergence=NA, loglik=NA) {
        super$initialize(model.par, variational.par, criteria, convergence, loglik)
        self$penalty <- penalty
      }
    )
)

## TODO: add function to for network representation
