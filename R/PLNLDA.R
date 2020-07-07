#' Poisson lognormal model towards Linear Discriminant Analysis
#'
#' Fit the Poisson lognormal for LDA with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param grouping a factor specifying the class of each observation used for discriminant analysis.
#' @param control a list for controlling the optimization process. See details.
#'
#' @return an R6 object with class [PLNLDAfit()]
#'
#' @details The parameter `control` is a list controlling the optimization with the following entries:
#' * "covariance" character setting the model for the covariance matrix. Either "full" or "spherical". Default is "full".
#' * "trace" integer for verbosity.
#' * "inception" Set up the initialization. By default, the model is initialized with a multivariate linear model applied on log-transformed data. However, the user can provide a PLNfit (typically obtained from a previous fit), which often speed up the inference.
#' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-8
#' * "ftol_abs" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0
#' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4
#' * "xtol_abs" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 0
#' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
#' * "maxtime" stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)
#' * "algorithm" the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS", "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".
#'
#'
#' @rdname PLNLDA
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
#' @seealso The class [`PLNLDAfit`]
#' @importFrom stats model.frame model.matrix model.response model.offset
#' @export
PLNLDA <- function(formula, data, subset, weights, grouping, control = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## look for grouping in the data or the parent frame
  if (class(try(eval(grouping), silent = TRUE)) == "try-error") {
    grouping <- try(eval(substitute(grouping), data), silent = TRUE)
    if (class(grouping) == "try-error") stop("invalid grouping")
  }
  grouping <- as.factor(grouping)

  ## treatment of the design, which is specific to LDA
  # - save the covariates
  covar <- args$X
  # - remove intercept so that 'grouping' describes group means
  xint <- match("(Intercept)", colnames(covar), nomatch = 0L)
  if (xint > 0L) covar <- covar[, -xint, drop = FALSE]
  # - build the design matrix encompassing covariates and LDA grouping
  design_group <- model.matrix(~ grouping + 0)
  args$X <- cbind(covar, design_group)

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- PLN_param(control, nrow(args$Y), ncol(args$Y), ncol(args$X))

  ## Initialize LDA by adjusting a PLN
  myLDA <- PLNLDAfit$new(grouping, args$Y, args$X, args$O, args$w,
                         args$model, args$xlevels, ctrl)

  ## Compute the group means
  if (ctrl$trace > 0) cat("\n Performing discriminant Analysis...")
  myLDA$optimize(args$X, covar, design_group, ctrl)

  ## Post-treatment: prepare LDA visualization
  myLDA$postTreatment(args$Y, args$X, args$O)

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myLDA
}
