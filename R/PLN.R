##' @title Fit a Poisson lognormal model
##'
##' @description two methods are available for specifying the models (with formulas or matrices)
##'
##' @param formula a formula
##' @param Y a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. Should include the intercept (a column of ones) if the default method is used.
##' @param O an optional (n x p) matrix of offsets.
##' @param control a list for controlling the optimization. See details.
##' @param Robject an R object, either a formula or a matrix
##' @param ... additional parameters. Not used
##'
##' @return an R6 object with class \code{\link[=PLNfit-class]{PLNfit}}
##'
##' @details The parameter \code{control} is a list controlling the optimization with the following entries
##' \itemize{
##'  \item{"xtol"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"ftol"}{stop when an optimization step changes the objective function by less than xtol multiplied by the absolute value of the parameter. Default is 1e-6}
##'  \item{"maxit"}{stop when the number of iteration exceeds maxiter. Default is 10000}
##'  \item{"lbvar"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-5.}
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
  ctrl <- list(ftol=1e-6, xtol=1e-4, maxit=10000, lbvar=1e-5, method="MMA", trace=1)
  ctrl[names(control)] <- control

  ## problem dimensions and constant
  n  <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  KY <-sum(.logfactorial(Y))

  ## ===========================================
  ## OPTIMIZATION
  ##
  if (ctrl$trace > 0) cat("\n Adjusting the standard PLN model.")

  ## Initialization
  ## glm-Poisson model for the regression parameters
  Theta <- do.call(rbind, lapply(1:p, function(j) coefficients(glm.fit(X, Y[, j], offset = O[,j], family = poisson()))))
  ## 0 mean and minimum variance for the variational parameters
  par0 <- c(Theta, rep(0, n*p), rep(10*ctrl$lbvar, n*p))
  lower.bound <- c(rep(-Inf, p*d), rep(-Inf, p*n), rep(ctrl$lbvar, n*p))

  ## Now optimize with NLOPTR
  opts <- list("algorithm"   = paste("NLOPT_LD",ctrl$method, sep="_"),
               "maxeval"     = ctrl$maxit,
               "xtol_rel"    = ctrl$xtol,
               "ftol_rel"    = ctrl$ftol,
               "print_level" = max(0,ctrl$trace-1))

  optim.out <- nloptr(par0, eval_f = fn_optim_PLN_profiled_Cpp, lb = lower.bound, opts = opts,
                      Y = Y, X = X, O = O, KY = KY)
  ## ===========================================
  ## POST-TREATMENT
  ##
  Theta <- matrix(optim.out$solution[1:(p*d)]          , p,d)
  M     <- matrix(optim.out$solution[p*(d)   + 1:(n*p)], n,p)
  S     <- matrix(optim.out$solution[p*(d+n) + 1:(n*p)], n,p)
  Sigma <- crossprod(M)/n + diag(colMeans(S))
  Omega <- solve(Sigma)

  rownames(Theta) <- colnames(Y); colnames(Theta) <- colnames(X)
  dimnames(S)     <- dimnames(Y)
  dimnames(M)     <- dimnames(Y)
  rownames(Omega) <- colnames(Omega) <- colnames(Y)
  rownames(Sigma) <- colnames(Sigma) <- colnames(Y)

  ## compute some criteria for evaluation
  J   <- - optim.out$objective
  BIC <- J - (p * d + p*(p+1)/2) * log(n)
  ICL <- BIC - .5*n*p *log(2*pi*exp(1)) - sum(log(S))

  return(PLNfit$new(model.par       = list(Omega = Omega, Sigma = Sigma, Theta = Theta),
                    variational.par = list(M = M, S = S),
                    criteria        = c(J = J, BIC = BIC, ICL = ICL),
                    convergence     = data.frame(status = optim.out$status,
                                                 objective = optim.out$objective,
                                                 iterations=optim.out$iterations)))
}
