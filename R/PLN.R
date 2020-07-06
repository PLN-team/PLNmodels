#' Poisson lognormal model
#'
#' Fit the multivariate Poisson lognormal model with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets, weights).
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which PLN is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param control a list for controlling the optimization. See details.
#'
#' @return an R6 object with class [`PLNfit`]
#'
#' @details The parameter `control` is a list controlling the optimization with the following entries:
#' * "covariance" character setting the model for the covariance matrix. Either "full", "diagonal" or "spherical". Default is "full".
#' * "trace" integer for verbosity.
#' * "inception" Set up the initialization. By default, the model is initialized with a multivariate linear model applied on log-transformed data, and with the same formula as the one provided by the user. However, the user can provide a PLNfit (typically obtained from a previous fit), which sometimes speeds up the inference.
#' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6 when n < p, 1e-8 otherwise.
#' * "ftol_abs" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0
#' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4
#' * "xtol_abs" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 0
#' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
#' * "maxtime" stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)
#' * "algorithm" the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS", "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".
#'
#'
#' @rdname PLN
#' @include PLNfit-class.R
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' @seealso The class [`PLNfit`]
#' @importFrom stats model.frame model.matrix model.response model.offset model.weights terms
#' @export
PLN <- function(formula, data, subset, weights, control = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and eventually overwrite them by user-defined parameters
  ctrl <- PLN_param(control, nrow(args$Y), ncol(args$Y), ncol(args$X))

  ## initialization
  if (ctrl$trace > 0) cat("\n Initialization...")
  myPLN <- PLNfit$new(args$Y, args$X, args$O, args$w,
                      args$model, args$xlevels, ctrl)

  ## optimization
  if (ctrl$trace > 0) cat("\n Adjusting a PLN model with", ctrl$covariance,"covariance model")
  myPLN$optimize(args$Y, args$X, args$O, args$w, ctrl)

  ## post-treatment
  if (ctrl$trace > 0) cat("\n Post-treatments...")
  myPLN$postTreatment(args$Y, args$X, args$O, args$w)

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myPLN
}
