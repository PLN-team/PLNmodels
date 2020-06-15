#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNMM}} produces a collection of models which are instances of object with class \code{PLNMMfit}.
#' A \code{PLNMMfit} (say, with k components) is itself a collection of k \code{PLNfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for ...
#'
#' @param cluster the number of clusters of the current model
#' @param components a list with cluster component element, each of whom is a \code{PLNfit}.
#' @param model_par a list with the matrices associated with the final estimated parameters of the mixture model: Theta (covariates), Sigma (latent covariance), mu (vector of means/centers) and pi (vector of cluster proportions)
#' @param posteriorProbabilities matrix of posterior probabilities of class belonging
#' @param mixtureParam vector of cluster proportions
#' @param loglik variational lower bound of the loglikelihood
#' @param BIC variational lower bound of the BIC
#' @param ICL variational lower bound of the ICL
#' @param R_squared approximated goodness-of-fit criterion
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNMM}}, the class \code{\link[=PLNMMfamily]{PLNMMfamily}}
PLNMMfit <-
  R6Class(classname = "PLNMMfit",
    private = list(
      comp       = NA, # list of mixture components (PLNfit)
      tau        = NA, # posterior probabilities of cluster belonging
      R2         = NA, # approximated goodness of fit criterion
      monitoring = NA  # a list with optimization monitoring quantities
    ),
    public  = list(
      initialize = function(inception, tau, monitoring=NA, R2=NA) {
        private$tau  <- tau
        components <- list()
        for (k_ in 1:ncol(tau)) components[[k_]] <- inception$clone()
        private$monitoring <- monitoring
        private$comp <- components
      },
      update = function(tau=NA, J=NA, R2=NA, monitoring=NA) {
        ## Theta=NA, Sigma=NA, M=NA, S=NA,
        ## later ... super$update(Theta, Sigma, M, S, J, R2, monitoring)
        if (!anyNA(tau))        private$tau    <- tau
        if (!anyNA(R2))         private$R2     <- R2
        if (!anyNA(monitoring)) private$monitoring <- monitoring
      },
      ## multi core ???
      optimize = function(responses, covariates, offsets, control) {

          ## ===========================================
          ## INITIALISATION
          tau  <- private$tau
          prop <- colMeans(private$tau)
          J_ic <- sapply(private$comp, function(comp) comp$loglik_vec)
          objective_old <- -sum(tau * J_ic) + sum(.xlogx(tau )) - self$n * sum(.xlogx(prop))

          ## ===========================================
          ## OPTIMISATION
          cond <- FALSE; iter <- 0
          objective   <- numeric(control$maxit_out)
          convergence <- numeric(control$maxit_out)
          while (!cond) {
            iter <- iter + 1
            if (control$trace > 1) cat("", iter)

            ## ---------------------------------------------------
            ## M - STEP
            ## UPDATE THE MIXTURE MODEL VIA OPTIMIZATION OF PLNMM
            for (k_ in seq.int(self$k)) {
              self$components[[k_]]$optimize(responses, covariates, offsets, tau[, k_], control)
            }
            prop <- colMeans(tau)

            ## ---------------------------------------------------
            ## E - STEP
            ## UPDATE THE POSTERIOR PROBABILITIES
            if (self$k > 1) { # only needed when at least 2 components!
              J_ic <- sapply(private$comp, function(comp) comp$loglik_vec)
              tau <- sweep(J_ic, 2, log(prop), "+") %>%
                apply(1, .softmax) %>% t() %>% .check_boundaries()
            }
            objective[iter]   <- -self$loglik
            convergence[iter] <- (objective[iter] - objective_old)/abs(objective[iter])

            ## Check convergence
            if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE
            objective_old <- objective[iter]

          }

          ## ===========================================
          ## OUTPUT
          ## formating parameters for output
          private$tau <- tau
          private$monitoring = list(objective        = objective[1:iter],
                                    convergence      = convergence[1:iter],
                                    outer_iterations = iter)
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Post treatment --------------------
      #' @description Update fields after optimization
      postTreatment = function(responses, covariates, offsets, weights, nullModel) {
        for (comp in self$components)
          comp$postTreatment(
            responses,
            covariates,
            offsets,
            weights,
            nullModel = nullModel
          )
      },
      show = function() {
        cat("Poisson Lognormal mixture model with",self$k,"components.\n")
        cat("* check fields $posteriorProb, $memberships, $mixtureParam and $components\n")
        cat("* each $component[[i]] is a PLNfit \n")
        cat("* average R2 is,",self$R_squared,"\n")
      },
      print = function() self$show()
    ),
    active = list(
      #' @field n number of samples
      n = function() {nrow(private$tau)},
      #' @field k number of components
      k = function() {ncol(private$tau)},
      #' @field components components of the mixture (PLNfits)
      components    = function() {private$comp},
      #' @field posteriorProb matrix ofposterior probability for cluster belonging
      posteriorProb = function() {private$tau},
      #' @field meberships vector for cluster index
      memberships   = function(value) {apply(private$tau, 1, which.max)},
      #' @field mixtureParam vector of cluster proporitions
      mixtureParam  = function() {colMeans(private$tau)},
      #' @field optim_par a list with parameters useful for monitoring the optimization
      optim_par  = function() {private$monitoring},
      #' @field nb_param number of parameters in the current PLN model
      nb_param      = function() {(self$k-1) + sum(sapply(self$components, function(model) model$nb_param))},
      #' @field entropy Entropy of the variational distribution of the cluster (multinomial)
      entropy       = function() {-sum(.xlogx(private$tau))},
      #' @field loglik variational lower bound of the loglikelihood
      loglik = function() {sum(self$loglik_vec)},
      #' @field loglik_vec element-wise variational lower bound of the loglikelihood
      loglik_vec = function() {
        J_ik <- sapply(private$comp, function(comp_) comp_$loglik_vec)
        rowSums(private$tau * J_ik) - rowSums(.xlogx(private$tau)) + rowSums(sweep(private$tau, 2, self$mixtureParam, "*"))
        },
      #' @field BIC variational lower bound of the BIC
      BIC        = function() {self$loglik - .5 * log(self$n) * self$nb_param},
      #' @field ICL variational lower bound of the ICL
      ICL        = function() {self$BIC - self$entropy},
      #' @field R_squared approximated goodness-of-fit criterion
      R_squared     = function() {sum(self$mixtureParam * sapply(self$components, function(model) model$R_squared))},
      #' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
      criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL, R_squared = self$R_squared)}

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
# PLNMMfit$set("public", "latent_pos",
# function(covariates, offsets) {
#   latentPos <- private$Mu + private$M + tcrossprod(covariates, private$Theta) + offsets
#   latentPos
# })

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

