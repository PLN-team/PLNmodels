##' Poisson lognormal model towards Principal Component Analysis
##'
##' Fit the PCA variants of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param control a list for controlling the optimization. See details.
##' @param ranks a vector of integer containing the successive ranks (or number of axes to be considered)
##' @param control.init a list for controling the optimization at initialization. See details.
##' @param control.main a list for controling the main optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNPCAfamily]{PLNPCAfamily}}, which contains
##' a collection of models with class \code{\link[=PLNPCAfit]{PLPCAfit}}
##'
##' @details The list of parameters \code{control.init} and \code{control.main} control the optimization of the intialization and the main process, with the following entries
##' \itemize{
##'   \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol_abs .}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol_rel multiplied by the absolute value of the parameter.}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol_abs.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "CCSAQ"  "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "TNEWTON_VAR1", "TNEWTON_VAR2". See NLOPTR documentation for further details. Default is "CCSAQ".}
##'  \item{"lbvar"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-5.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##'  \item{"inception"}{a PLNfit to start with. If NULL, a PLN is fitted on the . If an R6 object with class 'PLNfit' is given, it is used to initialize the model.}
##'  \item{"cores"}{The number of core used to paralellize jobs over the \code{ranks} vector. Default is 1.}
##' }
##'
##' @rdname PLNPCA
##' @examples
##' ## See the vignette: vignette("trichoptera", package="PLNmodels")
##' @seealso The classes \code{\link[=PLNnetworkfamily]{PLNPCAfamily}} and \code{\link[=PLNPCAfit]{PLNPCAfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNPCA <- function(formula, data, subset, ranks = 1:5,  control.init = list(), control.main = list()) {

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
  ctrl.init <- PLNPCA_param(control.init, n, p, "init")
  ctrl.main <- PLNPCA_param(control.main, n, p, "main")

  ## Instantiate the collection of PLN models, initialized by PLN with full rank
  if (ctrl.main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNPCAfamily$new(ranks = ranks, responses = Y, covariates = X, offsets = O, control = ctrl.init)

  ## Now adjust the PLN models
  if (ctrl.main$trace > 0) cat("\n\n Adjusting", length(ranks), "PLN models for PCA analysis.\n")
  myPLN$optimize(ctrl.main)

  ## Post-treatments: Compute pseudo-R2, rearrange criteria and the visualization for PCA
  if (ctrl.main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl.main$trace > 0) cat("\n DONE!\n")

  ## TODO formating the output to by (g)lm like
  ## TODO use the same output as in the broom package
  myPLN
}
