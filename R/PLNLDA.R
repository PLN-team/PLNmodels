##' Poisson lognormal model towards Linear Disciminant Analysis
##'
##' Fit the Poisson lognormal for LDA with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param grouping a factor specifying the class of each observation used for discriminant analysis.
##' @param control a list for controling the optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNLDAfit]{PLNLDAfit}}
##'
##' @details The list of parameters \code{control} tunes the optimization of the intialization and the main process, with the following entries
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
##' ## See the vignette: vignette("PLNLDA_trichoptera", package="PLNmodels")
##' @seealso The class \code{\link[=PLNLDAfit]{PLNLDAfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNLDA <- function(formula, data, subset, grouping, control = list(), ...) {

  ## create the call for the model frame
  call_frame <- match.call(expand.dots = FALSE)
  call_frame <- call_frame[c(1L, match(c("formula", "data", "subset"), names(call_frame), 0L))]
  call_frame[[1]] <- quote(stats::model.frame)

  ## eval the call in the parent environment
  frame <- eval(call_frame, parent.frame())

  ## create the set of matrices to fit the PLN model
  Y <- model.response(frame)
  n  <- nrow(Y); p <- ncol(Y) # problem dimension
  X <- model.matrix(terms(frame), frame)

  ## remove intercept in the covariates so that grouping describes means in each group
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L) X <- X[, -xint, drop = FALSE]
  d <- ncol(X)

  O <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- PLN_param(control, n, p)

  if (ctrl$trace > 0) cat("\n Initialization...")

  grouping <- as.factor(grouping)
  design_group <- model.matrix(~ grouping + 0)
  design_full <- cbind(X, design_group)

  myPLN <- PLN_internal(Y, design_full, O, ctrl)
  if (d > 0) {
    P <- (diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)) %*% myPLN$latent_pos(design_full, matrix(0, myPLN$n, myPLN$q))
    Group_Means <- t(rowsum(P, grouping) / tabulate(grouping))
  } else {
    Group_Means <- myPLN$model_par$Theta
  }
  colnames(Group_Means) <- colnames(design_group)

  if (ctrl$trace > 0) cat("\n Performing Discriminant Analysis...")
  myLDA <- PLNLDAfit$new(Theta = myPLN$model_par$Theta,
                         Group_Means = Group_Means,
                         Sigma = myPLN$model_par$Sigma,
                         grouping = grouping,
                         M = myPLN$var_par$M, S = myPLN$var_par$S,
                         J = myPLN$J, monitoring = myPLN$optim_par)
  myLDA$setVisualization()

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myLDA
}
