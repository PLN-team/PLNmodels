##' @title Fit a Poisson lognormal model
##'
##' @description two methods are available for specifing the models (with formulas or matrices)
##'
##' @param formula a formula
##' @param Y a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. SHould include the intercept (a column of one) if the default method is used.
##' @param O an optional (n x p) matrix of offsets.
##' @param control a list for controling the optimization. See details.
##' @param Robject an R object, either a formula or a matrix
##' @param ... additional parameters. Not used
##'
##' @return an R6 object with class \code{\link[=PLNfit-class]{PLNfit}}
##'
##' @details The parameter \code{control} is a list with the following entries
##' \itemize{
##'  \item{"factr"}{controls the L-BFGF-B procedure. See the documentation of \code{\link{optim}}. Default 1e8. Decrease if you experience instability or non monotonous J as a function of the rank}
##'  \item{"pgtol"}{controls the L-BFGF-B procedure. See the documentation of \code{\link{optim}}. Default 1e-2. Decrease if you experience instability or non monotonous J as a function of the rank}
##'  \item{"maxit"}{controls the L-BFGF-B procedure. See the documentaiton of \code{\link{optim}}. Default is 20000}
##'  \item{"lb.var"}{the minimum admissible value fr the variance parameter S in the variational approximation. Default is 1e-3.}
##'  \item{"cores"}{the number of cores. If Q has many entries, you might consider multiple cores. Default is 1.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##' }
##'
##' @rdname PLN
##' @include PLNfit-class.R
##' @examples
##' ## See the vignette TODO!!!
##' @seealso The class  \code{\link[=PLNfit-class]{PLNfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @importFrom nloptr nloptr
##' @export
PLN <- function(Robject, ...)
  UseMethod("PLN", Robject)

##' @rdname PLN
##' @export
PLN.formula <- function(formula, control = list()) {

  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLN.default(Y, X, O, control))
}

##' @rdname PLN
##' @export
PLN.default <- function(Y, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y)), control = list()) {

  ## ===========================================
  ## INITIALIZATION
  ##

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- list(ftol=1e-6, xtol=1e-4, maxit=10000, lbvar=1e-5, trace=1)
  ctrl[names(control)] <- control

  ## problem dimensions and constant
  n  <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  KY <-sum(.logfactorial(Y))

  fn_optim_PLN <- function(par) {
    Theta <- matrix(par[1:(p*d)]                         , p,d)
    M     <- matrix(par[p*d          + 1:(n*p)], n,p)
    S     <- matrix(par[(n+d)*p + 1:(n*p)], n,p)

    Omega <- n * chol2inv(chol(crossprod(M) + diag(colSums(S))))
    logDetOmega <- determinant(Omega, logarithm=TRUE)$modulus

    Z <- O + tcrossprod(X, Theta) + M
    A <- exp (.trunc(Z + .5*S))
    logP.Z  <- n/2 * (logDetOmega - sum(diag(Omega)*colMeans(S))) - .5*sum(diag(Omega %*% crossprod(M)))

    gr.Theta <- crossprod(X, A - Y)
    gr.M  <- M %*% Omega + A - Y
    gr.S  <- .5 * (matrix(rep(diag(Omega),n), n, p, byrow = TRUE) + A - 1/S)

    return(list(
      "objective" = sum(as.numeric(A - Y*Z)) - logP.Z - .5*sum(log(S)+1) + KY,
      "gradient"  = c(gr.Theta,gr.M,gr.S)
    ))
  }

  ## ===========================================
  ## OPTIMIZATION
  ##
  if (ctrl$trace > 0) cat("\n Adjusting the PLN model.")

  ## Initialization

  ## glm-Poisson model for the regression parameters
  Theta <- do.call(rbind, lapply(1:p, function(j) coefficients(glm.fit(X, Y[, j], offset = O[,j], family = poisson()))))
  ## 0 mean and minimum variance for the variational parameters
  par0 <- c(Theta, rep(0, n*p), rep(10*ctrl$lbvar, n*p))
  ## set box constraint for the variance parameters
  lower.bound <- c(rep(-Inf, p*(d+n)), rep(ctrl$lbvar, n*p))

  ## Now optimize with NLOPTR
  opts <- list("algorithm"   = "NLOPT_LD_MMA",
               "maxeval"     = ctrl$maxit,
               "xtol_rel"    = ctrl$xtol,
               "ftol_abs"    = ctrl$ftol,
               "ftol_rel"    = ctrl$ftol,
               "print_level" = ctrl$trace)

  optim.out <- nloptr(par0, eval_f = fn_optim_PLN, lb = lower.bound, opts = opts)

  ## ===========================================
  ## POST-TREATMENT
  ##
  Theta <- matrix(optim.out$solution[1:(p*d)]          ,d, p)
  M     <- matrix(optim.out$solution[p*d     + 1:(n*p)], n,p)
  S     <- matrix(optim.out$solution[(n+d)*p + 1:(n*p)], n,p)
  Sigma <- crossprod(M)/n + diag(colMeans(S))
  Omega <- solve(Sigma)
  colnames(Theta) <- colnames(Y); rownames(Theta) <- colnames(X)
  dimnames(S)     <- dimnames(Y)
  dimnames(M)     <- dimnames(Y)
  rownames(Omega) <- colnames(Omega) <- colnames(Y)

  ## compute some criteria for evaluation
  J   <- - optim.out$objective
  BIC <- J - (p * d + p*(p+1)/2) * log(n)
  ICL <- BIC - .5*n*p *log(2*pi*exp(1)) - sum(log(S))

  return(PLNfit$new(model.par       = list(Omega = Omega, Sigma = Sigma, Theta = Theta),
                      variational.par = list(M = M, S = S),
                      criteria        = c(J = J, BIC = BIC, ICL = ICL),
                      convergence     = data.frame(optim.out$message, optim.out$objective)))
}

