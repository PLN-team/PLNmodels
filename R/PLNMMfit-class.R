#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNMM}} produces a collection of models which are instances of object with class \code{PLNMMfit}.
#' A \code{PLNMMfit} (say, with k components) is itself a collection of k \code{PLNfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for ...
#'
#' @field cluster the number of clusters of the current model
#' @field components a list with cluster component element, each of whom is a \code{PLNfit}.
#' @field model_par a list with the matrices associated with the final estimated parameters of the mixture model: Theta (covariates), Sigma (latent covariance), mu (vector of means/centers) and pi (vector of cluster proportions)
#' @field posteriorProbabilities matrix of posterior probabilities of class belonging
#' @filed mixtureParam vector of cluster proportions
#' @field loglik variational lower bound of the loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field criteria a vector with loglik, BIC, ICL, R_squared and degrees of freedom
#' @field degrees_freedom number of parameters in the current model
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNMM}}, the class \code{\link[=PLNMMfamily]{PLNMMfamily}}
PLNMMfit <-
  R6Class(classname = "PLNMMfit",
    public  = list(
      initialize = function(inception, tau, J=NA, monitoring=NA, R2=NA) {
        private$tau  <- tau
        components <- list()
        for (k_ in 1:ncol(tau)) components[[k_]] <- inception$clone()
        private$J    <- J
        private$monitoring <- monitoring
        private$comp <- components
      },
      update = function(tau=NA, J=NA, R2=NA, monitoring=NA) {
        ## Theta=NA, Sigma=NA, M=NA, S=NA,
        ## later ... super$update(Theta, Sigma, M, S, J, R2, monitoring)
        if (!anyNA(tau))        private$tau    <- tau
        if (!anyNA(J))          private$J      <- J
        if (!anyNA(R2))         private$R2     <- R2
        if (!anyNA(monitoring)) private$monitoring <- monitoring
      }
    ),
    private = list(
      comp       = NA, # list of mixture components
      tau        = NA, # posterior probabilities of calss belonging
      R2         = NA, # approximated goodness of fit criterion
      J          = NA, # approximated loglikelihood
      monitoring = NA  # a list with optimization monitoring quantities
    ),
    active = list(
      n = function() {nrow(private$tau)},
      k = function() {ncol(private$tau)},
      components = function() {private$comp},
      posteriorProb = function() {private$tau},
      memberships = function(value) {apply(private$tau, 1, which.max)},
      mixtureParam = function() {colMeans(private$tau)},
      optim_par = function() {private$monitoring},
      degrees_freedom = function() {self$k + sum(sapply(self$components, function(model) model$degrees_freedom))},
      loglik    = function() {private$J},
      BIC       = function() {self$loglik - .5 * log(self$n) * self$degrees_freedom},
      ## ICL is wrong in this context
      ICL       = function() {self$BIC - .5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(private$S)))},
      R_squared = function() {private$R2},
      criteria  = function() {c(degrees_freedom = self$degrees_freedom, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL, R_squared = self$R_squared)}
    )
)

PLNMMfit$set("public", "optimize",
function(responses, covariates, offsets, tau, opts) {
  for (k_ in seq.int(self$k)) {
    optim_out <- optimization_PLN(
      c(self$components[[k_]]$model_par$Theta, self$components[[k_]]$var_par$M, self$components[[k_]]$var_par$S),
      responses,
      covariates,
      offsets,
      tau[, k_],
      opts
    )
    self$components[[k_]]$update(
      Theta = optim_out$Theta,
      Sigma = optim_out$Sigma,
      M     = optim_out$M,
      S     = optim_out$S,
      J     = -optim_out$objective,
      Ji    = optim_out$loglik,
      monitoring = optim_out[c("objective", "iterations", "status", "message")]
    )
  }



})

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
# PLNMMfit$set("public", "latent_pos",
# function(covariates, offsets) {
#   latentPos <- private$Mu + private$M + tcrossprod(covariates, private$Theta) + offsets
#   latentPos
# })

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

PLNMMfit$set("public", "show",
function() {
  super$show(paste0("Poisson Lognormal mxiture model with ",self$k," components.\n"))
  cat("* Additional fields for PLNMM\n")
  cat("    coming... \n")
  cat("* Additional methods for PLNMM\n")
  cat("    coming... \n")
})

