##' @title Fit a Poisson lognormal model towards LInear Disciminant Analysis
##'
##' @description two methods are available for specifing the models (with formulas or matrices)
##'
##' @param Robject either (n x p) matrix of count data (used with argument grouping) or a formula to describe the relationship between the count matrix and the covariates (apart from the grouping)
##' @param X an optional (n x d) matrix of covariates. A vector of intercept is included by default. Ignored when Robject is a formula.
##' @param O an optional (n x p) matrix of offsets. Ignored when Robject is a formula.
##' @param grouping a factor specifying the class fo< each observation used for disciminant analysis. Ignored when Robject is a formula.
##' @param control a list for controling the optimization process. See details.
##' @param ... additional parameters. Not used
##'
##' @return an R6 object with class \code{\link[=PLNLDAfit]{PLNLDAfit}}
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
##' @rdname PLNLDA
##' @examples
##' ## See the vignette: vignette("trichoptera", package="PLNmodels")
##' @seealso The class \code{\link[=PLNPCAfit]{PLNLDAfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNLDA <- function(Robject, ...)
  UseMethod("PLNLDA", Robject)

##' @rdname PLNLDA
##' @export
PLNLDA.formula <- function(Robject, grouping, control, ...) {

  formula <- Robject
  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)

  ## remove intercept in the covariates so that grouping describes means in each group
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L) X <- X[, -xint, drop = FALSE]
  O         <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNLDA.default(Y, grouping, X, O, control))
}

##' @rdname PLNLDA
##' @export
PLNLDA.default <- function(Robject, grouping, X = NULL, O = matrix(0, nrow(Y), ncol(Y)),  control = list(), ...) {

  Y <- Robject; rm(Robject) # no copy made
  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- PLN_param(control, nrow(Y), ncol(Y))
  Y <- as.matrix(Y)
  if (!is.null(X)) X <- as.matrix(X)

  ## stopifnot(nrow(Y) > ncol(Y))
  if (ctrl$trace > 0) cat("\n Initialization...")
  design_group <- model.matrix( ~as.factor(grouping)+0)
  design_full <- cbind(X, design_group)
  myPLN <- PLN(Y, design_full, O, control = ctrl)
  if (!is.null(X)) {
    P <- (diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)) %*% myPLN$latent_pos(design_full, matrix(0, myPLN$n, myPLN$q))
    Theta <- t(rowsum(P, grouping) / tabulate(grouping))
  } else {
    Theta <- myPLN$model_par$Theta
  }

  if (ctrl$trace > 0) cat("\n Performing Discriminant Analysis...")
  myLDA <- PLNLDAfit$new(Theta = Theta,
                         Sigma = myPLN$model_par$Sigma,
                         grouping = grouping,
                         M = myPLN$var_par$M, S = myPLN$var_par$S,
                         J = myPLN$J, monitoring = myPLN$monitoring)
  myLDA$setVisualization()

  if (ctrl$trace > 0) cat("\n DONE!\n")
  return(myLDA)
}
