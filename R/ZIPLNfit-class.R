#' An R6 Class to represent a ZIPLNfit
#'
#' @description The function [ZIPLN()] fit a model which is an instance of a object with class [`ZIPLNfit`].
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()],
#' [predict()], [vcov()] and [standard_error()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call and used for predictions.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param type approximation scheme to compute the fisher information matrix. Either `wald` (default) or `louis`. \code{type = "louis"} results in smaller confidence intervals.
#'
#' @inherit PLN details
#'
#' @rdname ZIPLNfit
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
ZIPLNfit <- R6Class(
  classname = "ZIPLNfit",
  inherit = PLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  public = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize a [`PLNfit`] model
    #' @importFrom stats lm.wfit lm.fit poisson residuals coefficients runif
    ## TODO: Once "set" is supported by Roxygen go back to external definition using
    ## PLNfit$set("public", "initialize", { ... })
    ## See https://github.com/r-lib/roxygen2/issues/931
    initialize = function(responses, covariates, offsets, weights, formula, xlevels, control) {
      super$initialize(responses, covariates, offsets, weights, formula, xlevels, control)
      delta <- 1* (responses == 0)
      private$zeros <- delta
      LOGREGs <- suppressWarnings(lapply(1:self$p, function(j) glm.fit(covariates, delta[, j], family = binomial(link = "logit")) ))
      private$Pi <- LOGREGs %>% map(fitted) %>% do.call(cbind, .)
      private$Theta0 <- do.call(rbind, lapply(LOGREGs, coefficients))
    },
    #' @description
    #' Update a [`ZIPLNfit`] object
    #' @param Theta  matrix of regression matrix (Poisson component)
    #' @param Theta0 matrix of regression matrix (Bernoulli component)
    #' @param Sigma  variance-covariance matrix of the latent variables
    #' @param M      matrix of mean vectors for the variational approximation
    #' @param S2     matrix of variance vectors for the variational approximation
    #' @param Ji     vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2     approximate R^2 goodness-of-fit criterion
    #' @param Z      matrix of latent vectors (includes covariates and offset effects)
    #' @param A      matrix of fitted values
    #' @param monitoring a list with optimization monitoring quantities
    #' @return Update the current [`PLNfit`] object
    update = function(Theta=NA, Theta0=NA, Sigma=NA, M=NA, S2=NA, Ji=NA, R2=NA, Z=NA, A=NA, monitoring=NA) {
      super$update(Theta=Theta, Sigma=Sigma, M=M, S2=S2, Ji=Ji, R2=R2, Z=Z, A=A, monitoring=monitoring)
      if (!anyNA(Theta0)) private$Theta0 <- Theta0
    },

   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ## Optimizers ----------------------------
   #' @description Call to the C++ optimizer and update of the relevant fields

   optimize = function(responses, covariates, offsets, weights, control) {

     args <- list(Y = responses, X = covariates, O = offsets, configuration = control)
     args$init_parameters <- list(Omega = NA, Pi = private$Pi, Theta0 = t(private$Theta0), Theta = t(private$Theta), M = private$M, S = sqrt(private$S2))
     optim_out <- do.call(optimize_zi, args)

     Z <- offsets + covariates %*% optim_out$parameters$Theta
     A <- exp(offsets + optim_out$parameters$M + .5 * optim_out$parameters$S^2)

     Ji <- optim_out$vloglik
     attr(Ji, "weights") <- weights

     self$update(
        Theta      = t(optim_out$parameters$Theta),
        Theta0     = t(optim_out$parameters$Theta0),
        Sigma      = solve(optim_out$parameters$Omega),
        M          = optim_out$parameters$M,
        S2         = (optim_out$parameters$S)**2,
        Z          = Z,
        A          = A,
        Ji         = Ji,
        monitoring = list(
          iterations = optim_out$nb_iter,
          message    = optim_out$stop_reason,
          objective  = optim_out$criterion)
      )

   },
   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Post treatment functions --------------
    #' @description Update R2 field after optimization
    set_R2 = function(responses, covariates, offsets, weights, nullModel = NULL) {
      private$R2 <- NA
    }
 ),
 private = list(
    Theta0 = NA, # the model parameters for the covariable ('0'/Bernoulli part)
    Pi     = NA, # probabilities for being observed
    zeros  = NA  # an indicator of the zeros
 ),
 active = list(
    #' @field model_par a list with the matrices of parameters found in the model (Theta, Sigma, plus some others depending on the variant)
    model_par  = function() {list(Theta = private$Theta, Theta0 = private$Theta0, Sigma = private$Sigma)},
    #' @field var_par a list with two matrices, M and S2, which are the estimated parameters in the variational approximation
    var_par    = function() {list(M = private$M, S2 = private$S2, Pi = private$Pi)},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$Pi * private$zeros + (1 - private$Pi) * private$A},
    #' @field nb_param number of parameters in the current ZIPLN model
    nb_param   = function() {super$nb_param + self$d * self$p}
 )
)
