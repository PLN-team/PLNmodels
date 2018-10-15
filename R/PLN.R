##' Poisson lognormal model
##'
##' Fit the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param control a list for controlling the optimization. See details.
##'
##' @return an R6 object with class \code{\link[=PLNfit]{PLNfit}}
##'
##' @details The parameter \code{control} is a list controlling the optimization with the following entries
##' \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "TNEWTON_VAR1", "TNEWTON_VAR2". See NLOPTR documentation for further details. Default is "MMA".}
##'  \item{"lbvar"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-5.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##'  \item{"inception"}{Set up the intialization. By default, the model is initialized with a multivariate linear model applied on log-transformed data. However, the user can provide a PLNfit (typically obtained from a previsous fit), which often speed up the inference.}
##' }
##'
##' @rdname PLN
##' @include PLNfit-class.R
##' @examples
##' ## See the vignette
##' @seealso The class  \code{\link[=PLNfit]{PLNfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset model.weights terms
##' @export
PLN <- function(formula, data, subset, weights, control = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and overwrite by user defined parameters
  args$ctrl <- PLN_param(control, nrow(args$Y), ncol(args$Y))

  ## call to the fitting function
  res <- do.call(PLN_internal, args)
  res
}

PLN_internal <- function(Y, X, O, w, ctrl) {

  ## ===========================================
  ## INITIALIZATION
  ##
  ## problem dimensions
  n  <- nrow(Y); p <- ncol(Y); d <- ncol(X)

  ## Check weights
  weighted <- ifelse(is.null(w), FALSE, TRUE)
  if (!weighted) w <- rep(1.0,n)

  ## get an initial point for optimization
  if (ctrl$trace > 0) cat("\n Initialization...")
  par0 <- initializePLN(Y, X, O, w, ctrl) # Theta, M, S

  ## ===========================================
  ## OPTIMIZATION
  ##
  par0 <- c(par0$Theta, par0$M, par0$S)

  if (ctrl$trace > 0) cat("\n Adjusting the standard PLN model")
  if (ctrl$trace > 0 & weighted) cat(" (weigthed)")
  opts <- list(
    "algorithm"   = ctrl$method,
    "maxeval"     = ctrl$maxeval,
    "ftol_rel"    = ctrl$ftol_rel,
    "ftol_abs"    = ctrl$ftol_abs,
    "xtol_rel"    = ctrl$xtol_rel,
    "xtol_abs"    = c(rep(0, p*d), rep(0, p*n), rep(ctrl$xtol_abs, n*p)),
    "lower_bound" = c(rep(-Inf, p*d), rep(-Inf, p*n), rep(ctrl$lbvar, n*p)),
    "weighted"    = weighted
  )

  optim_out <- optimization_PLN(par0, Y, X, O, w, opts)
  optim_out$message <- statusToMessage(optim_out$status)
  ## ===========================================
  ## POST-TREATMENT
  ##

  rownames(optim_out$Theta) <- colnames(Y)
  colnames(optim_out$Theta) <- colnames(X)
  rownames(optim_out$Sigma) <- colnames(optim_out$Sigma) <- colnames(Y)
  dimnames(optim_out$S) <- dimnames(optim_out$M) <- dimnames(Y)

  myPLN <- PLNfit$new(
    Theta = optim_out$Theta,
    Sigma = optim_out$Sigma,
    M     = optim_out$M,
    S     = optim_out$S,
    J     = -optim_out$objective,
    Ji    = optim_out$loglik,
    monitoring = optim_out[c("objective", "iterations", "status", "message")]
  )

  if (ctrl$trace > 0) cat("\n Computing (pseudo) R2")
  myPLN$set_R2(Y, X, O)

  if (ctrl$trace > 0) cat("\n DONE!\n")
  myPLN
}

## Extract the model used for initializing the whole family
#' @importFrom stats lm.wfit lm.fit poisson residuals coefficients runif
initializePLN <- function(Y, X, O, w, control) {

  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

  ## User defined (from a previous fit, for instance)
  if (isPLNfit(control$inception)) {
    if (control$trace > 1) cat("\n User defined inceptive PLN model")
    stopifnot(isTRUE(all.equal(dim(control$inception$model_par$Theta), c(p,d))),
              isTRUE(all.equal(dim(control$inception$var_par$M), c(n,p))),
              isTRUE(all.equal(dim(control$inception$var_par$S), c(n,p))))
    init <- list(Theta = control$inception$model_par$Theta,
                 M     = control$inception$var_par$M,
                 S     = control$inception$var_par$S)
  ## Default LM + log transformation
  } else {
    if (control$trace > 1) cat("\n Use LM after log transformation to define the inceptive model")
    LMs   <- lapply(1:p, function(j) lm.wfit(X, log(1 + Y[,j]), w, offset =  O[,j]) )
    Theta <- do.call(rbind, lapply(LMs, coefficients))
    M     <- do.call(cbind, lapply(LMs, residuals))
    S <- matrix(10 * control$lbvar, n, p)
    init <- list(Theta = Theta, M = M, S = S)
  }

  init
}

