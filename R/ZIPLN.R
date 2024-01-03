#' Zero Inflated Poisson lognormal model
#'
#' Fit the multivariate Zero Inflated Poisson lognormal model with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets, subset).
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. See details
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which PLN is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param lambda a positive number controlling the level of sparsity in the regression parameters of the PLN component. Default to 0 (no sparsity).
#' @param rho a positive number controlling the level of sparsity in the inverse covariance matrix. Default to 0 (no sparsity).
#' @param control a list for controlling the optimization. See details.
#'
#' @return an R6 object with class [`ZIPLNfit`]
#'
#' @details The parameter `control` is a list controlling the optimization with the following entries:
#' * covariance a character setting the model for the covariance matrix: either "full", "diagonal", "spherical" or "sparse". Default is "full", automatically set to sparse if rho > 0.
#' * "trace" integer for verbosity.
#' * "maxiter_out" control outer optimization (VEM) maximum number of iterations in the variational E-M algorithm. Default is 100.
#' * "ftol_out" control outer optimization (VEM) stop when an full V-EM iteration  changes the ELBO by less than ftol_out. Default is 1e-4.
#' * "ftol_rel" control inner optimization (VE step) stop when an optimization step changes in the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6 when n < p, 1e-8 otherwise.
#' * "ftol_abs" control inner optimization (VE step) stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0
#' * "xtol_rel" control inner optimization (VE step) stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4
#' * "xtol_abs" control inner optimization (VE step) stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 0
#' * "maxeval" control inner optimization (VE step) stop when the number of iteration exceeds maxeval. Default is 1000
#' * "maxtime" control inner optimization (VE step) stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)
#' * "algorithm" control inner optimization (VE step): the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS", "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".
#'
#' @rdname ZIPLN
#' @include ZIPLNfit-class.R
#' @examples
#' data(scRNA)
#' # data subsample: only 100 random cell and the 50 most varying transcript
#' scRNA        <- scRNA[sample.int(nrow(scRNA), 100), ]
#' scRNA$counts <- scRNA$counts[, 1:50]
#' myPLN_full   <- ZIPLN(counts ~ 1 + cell_line + offset(log(total_counts)) | 1, data = scRNA)
#' myPLN_sparse <- ZIPLN(counts ~ 1 + offset(log(total_counts)), rho = .5, data = scRNA)
#' myPLN_full$criteria    # better BIC with sparse version
#' myPLN_sparse$criteria
#' @seealso The class [`ZIPLNfit`]
#' @importFrom stats model.frame model.matrix model.response model.offset terms as.formula
#' @export
ZIPLN <- function(formula, data, subset, rho = 0, lambda = 0, control = list()) {

  ## extract the data matrices and weights
  args <- extract_model_zi(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and eventually overwrite them by user-defined parameters
  ctrl <- ZIPLN_param(control, nrow(args$Y), ncol(args$Y))
  ctrl$lambda  <- lambda
  ctrl$rho     <- rho
  ctrl$ziparam <- args$ziparam
  if (ctrl$rho > 0) ctrl$covariance <- "sparse"

  ## handling intercept term for penalized regression
  if (attr(terms(as.formula(args$formula)), "intercept") >0 & !ctrl$penalize_intercept)
    ctrl$ind_intercept <- match("(Intercept)", colnames(args$X))

  ## initialization
  if (ctrl$trace > 0) cat("\n Initialization...")
  myPLN <- ZIPLNfit$new(args$Y, args$X, args$O, args$formula, args$xlevels, ctrl)

  ## optimization
  if (ctrl$trace > 0) cat("\n Adjusting a ZI-PLN model with",
                          ctrl$covariance,"covariance model and",
                          ctrl$ziparam, "specific parameter(s) in Zero inflation component.")
  myPLN$optimize(args$Y, args$X, args$O, ctrl)

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myPLN
}

## -----------------------------------------------------------------
##  Series of setter to default parameters for user's main functions

available_algorithms <- c("MMA", "CCSAQ", "LBFGS", "VAR1", "VAR2", "TNEWTON", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART")

ZIPLN_param <- function(control, n, p) {
  xtol_abs    <- ifelse(is.null(control$xtol_abs)   , 0         , control$xtol_abs)
  ctrl <- list(
    "ftol_out"   = 1e-8,
    "maxit_out"  = 100 ,
    "algorithm"   = "CCSAQ",
    "maxeval"     = 1000  ,
    "maxtime"     = -1     ,
    "ftol_rel"    = 1e-8 ,
    "ftol_abs"    = 0,
    "xtol_rel"    = 1e-8,
    "xtol_abs"    = xtol_abs,
    "trace"       = 1,
    "covariance"  = "full",
    "penalize_intercept" = FALSE
  )
  ctrl[names(control)] <- control
  stopifnot(ctrl$algorithm %in% available_algorithms)
  ctrl
}

