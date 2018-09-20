##' Poisson lognormal model
##'
##' Fit the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
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
##' @importFrom stats model.frame model.matrix model.response model.offset terms
##' @export
PLN <- function(formula, data, subset, weights, control = list()) {

  ## create the call for the model frame
  call_frame <- match.call(expand.dots = FALSE)
  call_frame <- call_frame[c(1L, match(c("formula", "data", "subset", "weights"), names(call_frame), 0L))]
  call_frame[[1]] <- quote(stats::model.frame)

  ## eval the call in the parent environment
  frame <- eval(call_frame, parent.frame(2))

  ## create the set of matrices to fit the PLN model
  Y <- model.response(frame)
  X <- model.matrix(terms(frame), frame)
  O <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))
  w <- model.weights(frame)
  if (!is.null(w)) stopifnot(all(w > 0) && length(w) == nrow(Y))

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- PLN_param(control, nrow(Y), ncol(Y))

  ## call to the fitting function
  res <- PLN_internal(
    as.matrix(Y),
    as.matrix(X),
    as.matrix(O),
    as.vector(w),
    ctrl
    )

  res
}

PLN_internal <- function(Y, X, O, w, ctrl) {

  ## ===========================================
  ## INITIALIZATION
  ##
  ## problem dimensions
  n  <- nrow(Y); p <- ncol(Y); d <- ncol(X)

  ## get an initial point for optimization
  par0 <- initializePLN(Y, X, O, ctrl) # Theta, M, S

  ## ===========================================
  ## OPTIMIZATION
  ##
  if (ctrl$trace > 0) cat("\n Adjusting the standard PLN model.")

  par0 <- c(par0$Theta, par0$M, par0$S)

  opts <- list(
    "algorithm"   = ctrl$method,
    "maxeval"     = ctrl$maxeval,
    "ftol_rel"    = ctrl$ftol_rel,
    "ftol_abs"    = ctrl$ftol_abs,
    "xtol_rel"    = ctrl$xtol_rel,
    "xtol_abs"    = c(rep(0, p*d), rep(0, p*n), rep(ctrl$xtol_abs, n*p)),
    "lower_bound" = c(rep(-Inf, p*d), rep(-Inf, p*n), rep(ctrl$lbvar, n*p))
  )

  if (is.null(w)) {
    optim.out <- optimization_PLN(par0, Y, X, O, opts)
  } else {
    optim.out <- optimization_wPLN(par0, Y, X, O, w, opts)
  }

  optim.out$message <- statusToMessage(optim.out$status)

  ## ===========================================
  ## POST-TREATMENT
  ##
  Theta <- matrix(optim.out$solution[1:(p*d)]          , p,d)
  M     <- matrix(optim.out$solution[p*(d)   + 1:(n*p)], n,p)
  S     <- matrix(optim.out$solution[p*(d+n) + 1:(n*p)], n,p)
  Sigma <- crossprod(M)/n + diag(colMeans(S), nrow = p, ncol = p)

  rownames(Theta) <- colnames(Y); colnames(Theta) <- colnames(X)
  rownames(Sigma) <- colnames(Sigma) <- colnames(Y)
  dimnames(S)     <- dimnames(Y)
  dimnames(M)     <- dimnames(Y)

  myPLN <- PLNfit$new(
    Theta = Theta, Sigma = Sigma, M = M, S = S, J = -optim.out$objective,
    monitoring = optim.out[c("objective", "iterations", "status", "message")]
    )
  ## Compute R2
  loglik <- logLikPoisson(Y, myPLN$latent_pos(X, O))
  lmin   <- logLikPoisson(Y, nullModelPoisson(Y, X, O))
  lmax   <- logLikPoisson(Y, fullModelPoisson(Y))
  myPLN$update(R2 = (loglik - lmin) / (lmax - lmin))
  myPLN
}

## Extract the model used for initializing the whole family
#' @importFrom stats glm.fit lm.fit poisson residuals coefficients runif
initializePLN <- function(Y, X, O, control) {

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
    LMs  <- lapply(1:p, function(j) lm.fit(X, log(1 + Y[,j]), offset =  O[,j]) )
    Theta <- do.call(rbind, lapply(LMs, coefficients))
    M     <- do.call(cbind, lapply(LMs, residuals))
    S <- matrix(10 * control$lbvar, n, p)
    init <- list(Theta = Theta, M = M, S = S)
  }

  init
}

