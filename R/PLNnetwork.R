##' @title Fit a Poisson lognormal model towards network inference
##'
##' @description two methods are available for specifing the models (with formulas or matrices)
##'
##' @param formula a formula
##' @param Y a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. SHould include the intercept (a column of one) if the default method is used.
##' @param O an optional (n x p) matrix of offsets.
##' @param penalties a vector of positive real number controling the level of sparisty of the underlying network
##' @param control a list for controling the optimization. See details.
##' @param Robject an R object, either a formula or a matrix
##' @param ... additional parameters. Not used
##'
##' @return an R6 object with class \code{\link[=PLNnetworkfamily-class]{PLNnetworkfamily}}, which contains
##' a collection of models with class \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
##'
##' @details The parameter \code{control} is a list with the following entries
##' \itemize{
##'  \item{"out.tol"}{outer solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6}
##'  \item{"out.maxit"}{outer solver stops when the number of iteration exeeeds maxiter. Default is 10000}
##'  \item{"xtol"}{inner solver stops when an optimization step changes every parameters by less than xtol multiply by the absolute value of the parameter. Default is 1e-4}
##'  \item{"ftol"}{inner solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6}
##'  \item{"maxit"}{inner solver stops when the number of iteration exeeeds maxiter. Default is 10000}
##'  \item{"lbvar"}{the lower bound (box constraint) for the variational variance parameters. Default is .Machine$double.eps.}
##'  \item{"lbvar.init"}{the lower bound (box constraint) for the variational variance parameters for the unpenalized model. Default is 1e-5.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##' }
##'
##' @rdname PLNnetwork
##' @examples
##' ## See the vignette TODO
##' @seealso The classes \code{\link[=PLNnetworkfamily-class]{PLNnetworkfamily}} and \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNnetwork <- function(Robject, ...)
  UseMethod("PLNnetwork", Robject)

##' @rdname PLNnetwork
##' @export
PLNnetwork.formula <- function(formula, penalties = NULL,  control = list()) {

  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNnetwork.default(Y, X, O, penalties, control))
}

##' @rdname PLNnetwork
##' @export
PLNnetwork.default <- function(Y, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y)), penalties = NULL, control = list()) {

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- list(out.tol = 1e-5, out.maxit = 50, ftol=1e-6, xtol=1e-4, maxit=10000, nPenalties = 25, lbvar.init=1e-5, penalize.diagonal = FALSE, lbvar=.Machine$double.eps, trace=1)
  ctrl[names(control)] <- control
  ctrl.init <- list(ftol=ctrl$ftol, xtol=ctrl$xtol, maxit=ctrl$maxit, lbvar=ctrl$lbvar.init, trace=max(ctrl$trace,1))
  if (!is.null(penalties)) ctrl$nPenalties <- length(penalties)

  ## Instantiate the collection of PLN models
  if (ctrl$trace > 0) cat("\n Initialization")
  myPLN <- PLNnetworkfamily$new(nModels=ctrl$nPenalties, responses=Y, covariates=X, offsets=O, control=ctrl.init)

  ## Get an appropriate grid of penalties
  if (ctrl$trace > 0) cat("\n Recovering an appropriate grid of penalties.")
  myPLN$setPenalties(penalties, ctrl$nPenalties, ctrl$trace > 0)

  if (ctrl$trace > 0) cat("\n Adjusting", length(myPLN$penalties), "PLN models for sparse network inference.")
  myPLN$optimize(ctrl)
  if (ctrl$trace > 0) cat("\n DONE!\n")

  myPLN$setCriteria()

  return(myPLN)
}

