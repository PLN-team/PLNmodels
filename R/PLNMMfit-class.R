#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNMM}} produces a collection of models which are instances of object with class \code{PLNMMfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNMMfit_plot]{plot}}
#'
#' @field cluster the number of clusters of the current model
#' @field model_par a list with the matrices associated with the estimated parameters of the mixture model: Theta (covariates), Sigma (latent covariance), mu (vector of means/centers) and pi (vector of cluster proportions)
#' @field var_par a list with three matrices, M, S, and tau, which are the estimated parameters in the variational approximation
#' @field optim_par a list with parameters useful for monitoring the optimization
#' @field loglik variational lower bound of the loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field criteria a vector with loglik, BIC, ICL, R_squared and degrees of freedom
#' @field degrees_freedom number of parameters in the current PLN model
#' @field membership
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNMM}}, the class \code{\link[=PLNMMfamily]{PLNMMfamily}}
PLNMMfit <-
  R6Class(classname = "PLNMMfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(Theta=NA, Sigma=NA, mu=NA, pi=NA, M=NA, S=NA, tau=NA, J=NA, monitoring=NA) {
        super$initialize(Theta, Sigma, M, S, J, monitoring)
        private$k   <- ncol(tau)
        private$mu  <- mu
        private$pi  <- pi
        private$tau <- tau
      },
      update = function(Theta=NA, Sigma=NA, mu=NA, pi=NA, M=NA, S=NA, tau=NA, J=NA, R2=NA, monitoring=NA) {
        super$update(Theta, Sigma, M, S, J, R2, monitoring)
        if (!anyNA(mu))  private$mu  <- mu
        if (!anyNA(pi))  private$pi  <- pi
        if (!anyNA(tau)) private$tau <- tau
      }
    ),
    private = list(
      mu  = NULL,
      tau = NULL,
      pi  = NULL
    ),
    active = list(
      k = function() {ncol(private$tau)},
      posteriorProb = function() {private$tau},
      memberships = function(value) {apply(private$tau, 1, which.max)},
      mixtureParam = function() {private$pi},
      degrees_freedom = function() {self$p * (self$d + self$k)},
      model_par = function() {
        par <- super$model_par
        par$mu <- private$mu
        par$pi <- private$pi
        par
      },
       var_par = function() {
        par <- super$var_par
        par$tau <- private$tau
        par
      }
    )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR INTERNAL USE -> PLNfamily
## ----------------------------------------------------------------------
## Should only be accessed BY PLNfamily but R6 friend class don't exist

# Positions in the (euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#
# @name PLNfit_latent_pos
#
# @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
# @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
PLNMMfit$set("public", "latent_pos",
function(covariates, offsets) {
  latentPos <- private$mu + private$M + tcrossprod(covariates, private$Theta) + offsets
  latentPos
})

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

PLNMMfit$set("public", "show",
function() {
  super$show(paste0("Poisson Lognormal mxiture model with ",self$k,"components.\n"))
  cat("* Additional fields for PLNMM\n")
  cat("    coming... \n")
  cat("* Additional methods for PLNMM\n")
  cat("    coming... \n")
})

