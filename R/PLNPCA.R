##' @title Fit a Poisson lognormal model towards Principal Component Analysis
##'
##' @description two methods are available for specifing the models (with formulas or matrices)
##'
##' @param formula a formula
##' @param Y a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. Should include the intercept (a column of one) if the default method is used.
##' @param O an optional (n x p) matrix of offsets.
##' @param ranks a vector of integer containing the successive ranks (or number of axes to be considered)
##' @param control.init a list for controling the optimization. See details.
##' @param control.main a list for controling the optimization. See details.
##' @param Robject an R object, either a formula or a matrix
##' @param ... additional parameters. Not used
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
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "TNEWTON_VAR1", "TNEWTON_VAR2". See NLOPTR documentation for further details. Default is "MMA".}
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
PLNPCA <- function(Robject, ...)
  UseMethod("PLNPCA", Robject)

##' @rdname PLNPCA
##' @export
PLNPCA.formula <- function(formula, ranks = 1:5,  control.init = list(), control.main = list()) {

  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNPCA.default(Y, X, O, ranks, control.init, control.main))
}

##' @rdname PLNPCA
##' @export
PLNPCA.default <- function(Y, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y)), ranks = 1:5,  control.init = list(), control.main = list()) {

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl.init <- list(ftol_rel = 1e-6, ftol_abs = 1e-4, xtol_rel = 1e-4, xtol_abs = 1e-5, maxeval = 10000, method = "MMA", lbvar = 1e-4, trace = 1, inception = ifelse(ncol(Y) < 500, "PLN", "LM"))
  ctrl.main <- list(ftol_rel = 1e-8, ftol_abs = 1e-5, xtol_rel = 1e-4, xtol_abs = 1e-5, maxeval = 10000, method = "MMA", lbvar = 1e-5, trace = 1, cores = 1)

  ctrl.init[names(control.init)] <- control.init
  ctrl.main[names(control.main)] <- control.main

  ## Instantiate the collection of PLN models, initialized by PLN with full rank
  if (ctrl.main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNPCAfamily$new(ranks=ranks, responses=Y, covariates=X, offsets=O, control=ctrl.init)

  ## Now adjust the PLN models
  if (ctrl.main$trace > 0) cat("\n\n Adjusting", length(ranks), "PLN models for PCA analysis.")
  myPLN$optimize(ctrl.main)

  ## Post-treatments: Compute pseudo-R2, rearrange criteria and the visualization for PCA
  if (ctrl.main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl.main$trace > 0) cat("\n DONE!\n")
  return(myPLN)
}
