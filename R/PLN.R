##' @title Fit a Poisson lognormal model with a variational algorithm
##'
##' @description two methods are available for specifying the models (with formulas or matrices)
##'
##' @param Robject an R object, either a formula or a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. A vector of intercept is included by default. Ignored when Robject is a formula.
##' @param O an optional (n x p) matrix of offsets. Ignored when Robject is a formula.
##' @param control a list for controlling the optimization. See details.
##' @param ... additional parameters for S3 compatibility. Not used
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
##'  \item{"inception"}{How to setup the intialization: either with a PLNfit (typically obtained from a previsous fit), or a character string between "LM" and "GLM".
##'   If "LM", a log transformation is applied to Y then a linear model for initialization.
##'   If "GLM", a GLM Poisson is used. Default is "LM".}
##' }
##'
##' @rdname PLN
##' @include PLNfit-class.R
##' @examples
##' ## See the vignette
##' @seealso The class  \code{\link[=PLNfit]{PLNfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLN <- function(Robject, ...) UseMethod("PLN", Robject)

##' @rdname PLN
##' @export
PLN.formula <- function(Robject, control = list(), ...) {

  formula <- Robject
  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLN.default(Y, X, O, control))
}

##' @rdname PLN
##' @export
PLN.default <- function(Robject, X = NULL, O = NULL, control = list(), ...) {

  Y <- as.matrix(Robject); rm(Robject)
  ## ===========================================
  ## INITIALIZATION
  ##
  ## problem dimensions
  n  <- nrow(Y); p <- ncol(Y)
  if (is.null(X)) X <- matrix(1, n, 1)
  if (is.null(O)) O <- matrix(0, n, p)
  d <- ncol(X)

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- PLN_param(control, n, p)

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

  optim.out <- optimization_PLN(par0, Y, X, O, opts)
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
  if(isTRUE(all.equal(tail(class(control$inception), 2), c("PLNfit", "R6")))) {
    if (control$trace > 1) cat("\n User defined inceptive PLN model")
    stopifnot(isTRUE(all.equal(dim(control$inception$model_par$Theta), c(p,d))),
              isTRUE(all.equal(dim(control$inception$var_par$M), c(n,p))),
              isTRUE(all.equal(dim(control$inception$var_par$S), c(n,p))))
    return(list(Theta = control$inception$model_par$Theta,
                M     = control$inception$var_par$M,
                S     = control$inception$var_par$S))

    ## GLM Poisson
  } else if (isTRUE(all.equal(is.character(control$inception), control$inception == "GLM"))) {
    if (control$trace > 1) cat("\n Use GLM Poisson to define the inceptive model")
    GLMs  <- lapply(1:p, function(j) glm.fit(X, Y[, j], offset = O[,j], family = poisson()))
    Theta <- do.call(rbind, lapply(GLMs, coefficients))
    M     <- do.call(cbind, lapply(GLMs, residuals, "deviance"))
    S <- matrix(10 * control$lbvar, n, p)
    return(list(Theta = Theta, M = M, S = S))

    ## default LM + log transformation
  } else {
    if (control$trace > 1) cat("\n Use LM after log transformation to define the inceptive model")
    LMs  <- lapply(1:p, function(j) lm.fit(X, log(1 + Y[,j]), offset =  O[,j]) )
    Theta <- do.call(rbind, lapply(LMs, coefficients))
    M     <- do.call(cbind, lapply(LMs, residuals))
    S <- matrix(10 * control$lbvar, n, p)
    return(list(Theta = Theta, M = M, S = S))
  }

}

