##' Poisson lognormal model towards Linear Disciminant Analysis
##'
##' Fit the Poisson lognormal for LDA with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param grouping a factor specifying the class of each observation used for discriminant analysis.
##' @param control a list for controling the optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNLDAfit]{PLNLDAfit}}
##'
##' @details The list of parameters \code{control} tunes the optimization of the intialization and the main process, with the following entries
##' \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
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
PLNLDA <- function(formula, data, subset, weights, grouping, control = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## treatment of the design, which is specific to LDA
  # - save the covariates
  covar <- args$X
  # - remove intercept so that 'grouping' describes group means
  xint <- match("(Intercept)", colnames(covar), nomatch = 0L)
  if (xint > 0L) covar <- covar[, -xint, drop = FALSE]
  # - build the design matrix encompassing covariates and LDA grouping
  grouping <- as.factor(grouping)
  design_group <- model.matrix(~ grouping + 0)
  args$X <- cbind(covar, design_group)

  ## define default control parameters for optim and overwrite by user defined parameters
  args$ctrl <- PLN_param(control, nrow(args$Y), ncol(args$Y), ncol(args$X))

  ## call to the fitting function
  myPLN <- do.call(PLN_internal, args)

  ## extract group means
  if (ncol(covar) > 0) {
    proj_orth_X <- (diag(myPLN$n) - covar %*% solve(crossprod(covar)) %*% t(covar))
    P <- proj_orth_X %*% myPLN$latent_pos(args$X, matrix(0, myPLN$n, myPLN$q))
    Mu <- t(rowsum(P, grouping) / tabulate(grouping))
  } else {
    Mu <- myPLN$model_par$Theta
  }
  colnames(Mu) <- colnames(design_group)

  ## vizualization for discriminant analysis
  if (args$ctrl$trace > 0) cat("\n Performing discriminant Analysis...")
  myLDA <- PLNLDAfit$new(
    Theta      = myPLN$model_par$Theta,
    Mu         = Mu,
    Sigma      = myPLN$model_par$Sigma,
    grouping   = grouping,
    M          = myPLN$var_par$M,
    S          = myPLN$var_par$S,
    J          = myPLN$J,
    monitoring = myPLN$optim_par
  )
  myLDA$setVisualization()

  if (args$ctrl$trace > 0) cat("\n DONE!\n")
  myLDA
}
