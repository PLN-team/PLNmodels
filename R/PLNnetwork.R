##' Poisson lognormal model towards sparse network inference
##'
##' Fit the sparse inverse covariance variant of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
##' @param penalties an optional vector of positive real number controling the level of sparisty of the underlying network. if NULL (the default), will be set internally
##' @param control_init a list for controling the optimization of the PLN model used at initialization. See details.
##' @param control_main a list for controling the main optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}}, which contains
##' a collection of models with class \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
##'
##' @details The list of parameters \code{control_main} controls the optimization of the main process, with the following entries
##'  \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol_abs .}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol_rel multiplied by the absolute value of the parameter.}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol_abs.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"algorithm"}{the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS",
##'     "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".}
##'  \item{"cores"}{integer for number of cores used. Default is 1.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##'  \item{"ftol_out"}{outer solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6}
##'  \item{"maxit_out"}{outer solver stops when the number of iteration exceeds out.maxit. Default is 50}
##'  \item{"penalize_diagonal"}{boolean: should the diagonal terms be penalized in the graphical-Lasso? Default is FALSE.}
##'  \item{"penalty_weights"}{p x p matrix of weights (defaut filled with 1) to adapt the amount of shrinkage to each pairs of node. Must be symmetric with positive values.}
##' }
##'
##' The list of parameters \code{control_init} controls the optimization process in the initialization and in the function \code{\link[=PLN]{PLN}}, plus two additional parameters:
##' \itemize{
##'  \item{"nPenalties"}{an integer that specified the number of values for the penalty grid when internally generated. Ignored when penalties is non NULL}
##'  \item{"min.ratio"}{the penalty grid ranges from the minimal value that produces a sparse to this value multiplied by \code{min.ratio}. Default is 0.1.}
##' }
##'
##' @rdname PLNnetwork
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
##' @seealso The classes \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}} and \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNnetwork <- function(formula, data, subset, weights, penalties = NULL, control_init = list(), control_main = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl_main <- PLNnetwork_param(control_main, nrow(args$Y), ncol(args$Y),ncol(args$X), weighted = !missing(weights))
  if (is.null(control_init$trace)) control_init$trace <- 0
  ctrl_init <- PLN_param(control_init, nrow(args$Y), ncol(args$Y), ncol(args$X))
  if (is.null(ctrl_init$nPenalties)) ctrl_init$nPenalties <- 30
  if (is.null(ctrl_init$min.ratio)) ctrl_init$min.ratio   <- .1
  ctrl_init$penalty_weights <- ctrl_main$penalty_weights

  ## Instantiate the collection of models
  if (ctrl_main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNnetworkfamily$new(penalties, args$Y, args$X, args$O, args$w, args$model, ctrl_init)

  ## Optimization
  if (ctrl_main$trace > 0) cat("\n Adjusting", length(myPLN$penalties), "PLN with sparse inverse covariance estimation\n")
  if (ctrl_main$trace) cat("\tJoint optimization alternating gradient descent and graphical-lasso\n")
  myPLN$optimize(ctrl_main)

  ## Post-treatments
  if (ctrl_main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl_main$trace > 0) cat("\n DONE!\n")
  myPLN
}
