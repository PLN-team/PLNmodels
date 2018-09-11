##' Poisson lognormal model towards sparse network inference
##'
##' Fit the sparse inverse covariance variant of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param control a list for controlling the optimization. See details.
##' @param penalties an optional vector of positive real number controling the level of sparisty of the underlying network. if NULL (the default), will be set internally
##' @param control.init a list for controling the optimization at initialization. See details.
##' @param control.main a list for controling the main optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}}, which contains
##' a collection of models with class \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
##'
##' @details The list of parameters \code{control.init} and \code{control.main} control the optimization of the initialization and the main process.
##'
##'  The following entries are shared by both \code{control.init} and \code{control.main} and mainly concern the optimization parameters of NLOPT. There values can be different in \code{control.init} and \code{control.main}
##'  \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol_abs .}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol_rel multiplied by the absolute value of the parameter.}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol_abs.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "VAR1", "VAR2". See NLOPTR documentation for further details. Default is "MMA".}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##' }
##'
##' The following parameter are specific to initialization
##' \itemize{
##'  \item{"inception"}{a PLNfit to start with. If NULL, a PLN is fitted. If an R6 object with class 'PLNfit' is given, it is used to initialize the model.}
##'  \item{"nPenalties"}{an integer that specified the number of values for the penalty grid when internally generated. Ignored when penalties is non NULL}
##'  \item{"min.ratio"}{the penalty grid ranges from the minimal value that produces a sparse to this value multiplied by \code{min.ratio}. Default is 0.01 for high dimensional problem, 0.001 otherwise.}
##' }
##'
##' The following parameter are specific to main iterative process
##' \itemize{
##'  \item{"ftol_out"}{outer solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6}
##'  \item{"maxit_out"}{outer solver stops when the number of iteration exceeds out.maxit. Default is 50}
##'  \item{"penalize.diagonal"}{boolean: should the diagonal terms be penalized in the graphical-Lasso? Default is FALSE.}
##' }
##'
##' @rdname PLNnetwork
##' @examples
##' ## See the vignette
##' @seealso The classes \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}} and \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNnetwork <- function(formula, data, subset, penalties = NULL, control.init = list(), control.main = list()) {

  ## create the call for the model frame
  call_frame <- match.call(expand.dots = FALSE)
  call_frame <- call_frame[c(1L, match(c("formula", "data", "subset"), names(call_frame), 0L))]
  call_frame[[1]] <- quote(stats::model.frame)

  ## eval the call in the parent environment
  frame <- eval(call_frame, parent.frame())

  ## create the set of matrices to fit the PLN model
  Y <- model.response(frame)
  n  <- nrow(Y); p <- ncol(Y) # problem dimensions
  X <- model.matrix(terms(frame), frame)
  O <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl.init <- PLNnetwork_param(control.init, n, p, "init")
  ctrl.main <- PLNnetwork_param(control.main, n, p, "main")

  ## Instantiate the collection of PLN models
  if (ctrl.main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNnetworkfamily$new(penalties = penalties, responses = Y, covariates = X, offsets = O, control = ctrl.init)

  ## Main optimization
  if (ctrl.main$trace > 0) cat("\n Adjusting", length(myPLN$penalties), "PLN with sparse inverse covariance estimation\n")
  if (ctrl.main$trace) cat("\tJoint optimization alternating gradient descent and graphical-lasso\n")
  myPLN$optimize(ctrl.main)

  ## Post-treatments: compute pseudo-R2
  if (ctrl.main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl.main$trace > 0) cat("\n DONE!\n")
  myPLN
}
