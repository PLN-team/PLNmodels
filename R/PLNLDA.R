##' Poisson lognormal model towards Linear Disciminant Analysis
##'
##' Fit the Poisson lognormal for LDA with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
##' @param grouping a factor specifying the class of each observation used for discriminant analysis.
##' @param control a list for controling the optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNLDAfit]{PLNLDAfit}}
##'
##' @details The parameter \code{control} is a list controlling the optimization with the following entries
##' \itemize{
##'  \item{"covariance"}{character setting the model for the covariance matrix. Either "full" or "spherical". Default is "full".}
##'  \item{"trace"}{integer for verbosity.}
##'  \item{"inception"}{Set up the intialization. By default, the model is initialized with a multivariate linear model applied on log-transformed data. However, the user can provide a PLNfit (typically obtained from a previsous fit), which often speed up the inference.}
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"maxtime"}{stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)}
##'  \item{"algorithm"}{the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS",
##'     "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".}
##'  \item{"lower_bound"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-4.}
##' }
##'
##' @rdname PLNLDA
##' @examples
##' data(trichoptera)
##' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
##' myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
##' @seealso The class \code{\link[=PLNLDAfit]{PLNLDAfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
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
  myLDA <- PLNLDAfit$new(grouping, args$Y, args$X, args$O, args$w, args$model, ctrl)

  ## Compute the group means
  if (ctrl$trace > 0) cat("\n Performing discriminant Analysis...")
  myLDA$optimize(args$X, covar, design_group, ctrl)

  ## Post-treatment: prepare LDA vizualization
  myLDA$postTreatment(args$Y, args$X, args$O)

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myLDA
}
