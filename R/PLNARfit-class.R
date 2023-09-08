#' An R6 Class to represent a PLNARfit in a standard, general framework
#'
#' @description The function [PLN()] fit a model which is an instance of a object with class [`PLNARfit`].
#' Objects produced by the functions [PLNnetwork()], [PLNPCA()], [PLNmixture()] and [PLNLDA()] also enjoy the methods of [PLNARfit()] by inheritance.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()], [predict()], [vcov()] and [standard_error()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNARfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which PLN is called.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list-like structure for controlling the fit, see [PLN_param()].
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
#' @inherit PLN details
#'
#' @rdname PLNARfit
#' @include PLNARfit-class.R
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNARfit <- R6Class(
  classname = "PLNARfit",
  inherit = PLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## PRIVATE INTERNAL FIELDS
    Phi = NA, # the autoregressive matrix parameter

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torch_omega_squared_norm = function(x, Omega) {
      torch_sum(torch_mm(x, Omega) * x, dim = 2)
    },

    torch_frobenius = function(A, B) {
      torch_sum(A * B, dim = 2)
    },

    torch_elbo = function(data, params, index = torch_tensor(1:self$n)) {
      S2 <- torch_square(params$S[index])
      M_shift <- torch_roll(params$M)
      M_diff - params$M - torch_mm(params$Phi, M_shift)
      mu <- data$O[index] + torch_mm(data$X[index], params$B)
      mu_eps <- mu - torch_mm(Phi, mu)
      Omega_eps <- torch_inverse(params$Sigma - torch_mm(torch_mm(params$Phi, params$Sigma), torch_transpose(params$Phi)))
      
      res <- torch_exp(params$M[index] + .5 * S2) -
        sum(data$Y[index] * params$M[index]) -
        .5 * torch_logdet(params$Omega[index]) +
        .5 * torch_omega_squared_norm(torch_index_select(M, 2, 1:self$p) - mu, params$Omega) +
        .5 * torch_frobenius(torch_index_select(S2, 2, 1:self$p), torch_diag(params$Omega), dim = 2) -
        .5 * sum(data$w[index] - 1) * torch_logdet(params$Omega_eps[index]) +
        .5 * torch_omega_squared_norm(M_diff - mu_eps, Omega_eps) +
        .5 * torch_frobenius(torch_index_select(S2, 2, 2:data$w[index]), torch_diag(Omega_eps)) -
        .5 * torch_log(S2)
      res
    },

    torch_vloglik = function(data, params) {
      S2    <- torch_square(params$S)
      Ji <- - torch_elbo(data, params) + rowSums(.logfactorial(as.matrix(data$Y))) - .5 * self$p  
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF PRIVATE METHODS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize a [`PLNARfit`] object
    initialize = function(grouping, responses, covariates, offsets, weights, formula, control) {
      private$grouping <- grouping
      super$initialize(responses, cbind(covariates, model.matrix(~ grouping + 0)), offsets, weights, formula, control)
    }
  )
)
