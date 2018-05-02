##' @title Fit a Poisson lognormal model towards network inference
##'
##' @description two methods are available for specifying the models (with formulas or matrices)
##'
##' @param Robject an R object, either a formula or a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. A vector of intercept is included by default. Ignored when Robject is a formula.
##' @param O an optional (n x p) matrix of offsets. Ignored when Robject is a formula.
##' @param penalties an optional vector of positive real number controling the level of sparisty of the underlying network. if NULL (the default), will be set internally
##' @param control.init a list for controling the optimization at initialization. See details.
##' @param control.main a list for controling the main optimization process. See details.
##' @param ... additional parameters for S3 compatibility. Not used
##' @param ... additional parameters. Not used
##'
##' @return an R6 object with class \code{\link[=PLNnetworkfamily-class]{PLNnetworkfamily}}, which contains
##' a collection of models with class \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
##'
##' @details The list of parameters \code{control.init} and \code{control.main} control the optimization of the initialization and the main process.
##'
##'  The following entries are shared by both \code{control.init} and \code{control.main} and mainly concern the optimization parameters of NLOPT. There values can be different in \code{control.init} and \code{control.main}
##'  \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol_abs .}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol_rel multiplied by the absolute value of the parameter.}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol_abs.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "VAR1", "VAR2". See NLOPTR documentation for further details. Default is "MMA".}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##' }
##'
##' The following parameter are specific to initialization
##' \itemize{
##'  \item{"inception"}{a PLNfit to start with. If NULL, a PLN is fitted. If an R6 object with class 'PLNfit' is given, it is used to initialize the model.}
##'  \item{"nPenalties"}{an integer that specified the number of values for the penalty grid when internally generated. Ignored when penalties is non NULL}
##'  \item{"min.ratio"}{the penalty grid ranges from the minimal value that produces a sparse to this value multiplied by \code{min.ratio}. Default is 0.01 for high dimensional problem, 0.001 otherwise.}
##' }
##'
##' The following parameter are specific to main iterative process
##' \itemize{
##'  \item{"ftol_out"}{outer solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6}
##'  \item{"maxit_out"}{outer solver stops when the number of iteration exceeds out.maxit. Default is 50}
##'  \item{"penalize.diagonal"}{boolean: should the diagonal terms be penalized in the graphical-Lasso? Default is FALSE.}
##' }
##'
##' @rdname PLNnetwork
##' @examples
##' ## See the vignette
##' @seealso The classes \code{\link[=PLNnetworkfamily-class]{PLNnetworkfamily}} and \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNnetwork <- function(Robject, ...) UseMethod("PLNnetwork", Robject)

##' @rdname PLNnetwork
##' @export
PLNnetwork.formula <- function(Robject, penalties = NULL, approx = FALSE, control.init = list(), control.main = list(), ...) {

  formula <- Robject
  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNnetwork.default(Y, X, O, penalties, approx, control.init, control.main))
}

##' @rdname PLNnetwork
##' @export
PLNnetwork.default <- function(Robject, X = matrix(1, nrow = nrow(Robject)), O = matrix(0, nrow(Robject), ncol(Robject)),
                               penalties = NULL, approx=FALSE, control.init = list(), control.main=list(), ...) {

  Y <- Robject; rm(Robject) # no copy made
  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl.init <- PLNnetwork_param(control.init, nrow(Y), ncol(Y), "init")
  ctrl.main <- PLNnetwork_param(control.main, nrow(Y), ncol(Y), "main")

  ## approximation can be obtained by performing just one iteration in the joint optimization algorithm
  if (approx) ctrl.main$maxit_out <- 1

  ## Instantiate the collection of PLN models
  if (ctrl.main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNnetworkfamily$new(penalties = penalties, responses = Y, covariates = X, offsets = O, control = ctrl.init)

  ## Main optimization
  if (ctrl.main$trace > 0) cat("\n Adjusting", length(myPLN$penalties), "PLN with sparse inverse covariance estimation\n")
  if (ctrl.main$trace & approx) cat("\tTwo-step approach applying Graphical-Lasso on the inceptive PLN fit.\n")
  if (ctrl.main$trace & !approx) cat("\tJoint optimization alternating gradient descent and graphical-lasso\n")
  myPLN$optimize(ctrl.main)

  ## Post-treatments: compute pseudo-R2
  if (ctrl.main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl.main$trace > 0) cat("\n DONE!\n")
  return(myPLN)
}

##' @title Performs stability selection for Poisson lognormal models towards network inference
##'
##' @description two methods are available for specifying the models (with formulas or matrices)
##'
##' @param Robject an R object, either a formula or a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. A vector of intercept is included by default. Ignored when Robject is a formula.
##' @param O an optional (n x p) matrix of offsets. Ignored when Robject is a formula.
##' @param penalties an optional vector of positive real number controling the level of sparisty of the underlying network. if NULL (the default), will be set internally
##' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines th number of subsamples used in the stability selection. Automatically set to 100 subsamples with size \code{10*sqrt(n)} if \code{n >= 144} and \code{0.8*n} otherwise following Liu et al. (2010) recommandations.
##' @param approx a boolean for the type of optimization. if \code{FALSE}, perform the full alternating optimization scheme. if \code{TRUE} the fastest (yet approximated) two-step approach is used, first estimating a PLN model then applying graphical-Lasso on a grid of penalties. Default to TRUE.
##' @param control.init a list for controling the optimization of the initialization, that fits a standard PLN model with the \code{\link[=PLN]{PLN}} function. See details.
##' @param control.main a list for controling the optimization. See details.
##' @param mc.cores the number of cores to used. Default is 1.
##' @param ... additional parameters. Not used
##'
##' @return a list with the vector of penalty used, a p*(p-1)/2 x length(penalties) matrix of estimated probabilities of selection of the edges, and the list of subsamples
##'
##' @details See \code{\link[=PLNnetwork]{PLNnetwork}} for details about \code{control.main} and \code{control.init}
##'
##' @export
PLNnetwork_stabs <- function(Robject, ...) UseMethod("PLNnetwork_stabs", Robject)

##' @rdname PLNnetwork_stabs
##' @export
PLNnetwork_stabs.formula <- function(Robject, penalties = NULL, subsamples = NULL, approx = TRUE, control.init = list(), control.main = list(), mc.cores = 1, ...) {

  formula <- Robject
  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNnetwork_stabs.default(Y, X, O, penalties, subsamples, approx, control.init, control.main, mc.cores))
}

##' @rdname PLNnetwork_stabs
##' @export
PLNnetwork_stabs.default <- function(Robject, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y)),
                               penalties = NULL, subsamples = NULL, approx = TRUE, control.init = list(nPenalties = 10), control.main = list(), mc.cores = 1, ...) {

  Y <- Robject; rm(Robject) # no copy made
  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl.init <- PLNnetwork_param(control.init, nrow(Y), ncol(Y), "init")
  ctrl.init$trace <- 0
  ctrl.main <- PLNnetwork_param(control.main, nrow(Y), ncol(Y), "main")
  ctrl.main$trace <- 0

  ## get the common inceptive model to save time
  if ("PLNfit" %in% class(ctrl.init$inception)) {
    inception <- ctrl.init$inception
  } else {
    inception <- PLN(Y, X, O)
  }

  ## approximation can be obtained by performing just one iteration in the joint optimization algorithm
  if (approx) ctrl.main$maxit_out <- 1

  ## select default subsamples according
  if (is.null(subsamples)) {
    n <- nrow(Y)
    subsample.size <- round(ifelse(n >= 144, 10*sqrt(n), 0.8*n))
    subsamples <- replicate(20, sample.int(nrow(Y), subsample.size), simplify = FALSE)
  }

  ## Get an appropriate grid of penalties
  if (is.null(penalties)) {
    if (ctrl.init$trace > 1) cat("\n Recovering an appropriate grid of penalties.")
    max_pen <- max(abs(inception$model_par$Sigma[upper.tri(inception$model_par$Sigma)]))
    penalties <- 10^seq(log10(max_pen), log10(max_pen*ctrl.init$min.ratio), len = ctrl.init$nPenalties)
  } else {
    if (ctrl.init$trace > 1) cat("\nPenalties already set by the user")
    stopifnot(all(penalties > 0))
  }


  ## got for stability selection
  cat("\nStability Selection for PLNnetwork: ")
  cat("\nsubsampling: ")
  stabs_out <- mclapply(subsamples, function(subsample) {
    cat("+")
    inception_ <- inception$clone()
    inception_$update(M = inception$var_par$M[subsample, , drop = FALSE])
    inception_$update(S = inception$var_par$S[subsample, , drop = FALSE])
    ctrl.init$inception <- inception_
    myPLN <- PLNnetworkfamily$new(penalties  = penalties,
                                  responses  = Y[subsample, , drop = FALSE ],
                                  covariates = X[subsample, , drop = FALSE],
                                  offsets    = O[subsample, , drop = FALSE ], control = ctrl.init)
    myPLN$optimize(ctrl.main)
    nets <- do.call(cbind, lapply(myPLN$models, function(model) {
      as.matrix(model$latent_network())[upper.tri(diag(ncol(Y)))]
    }))
    nets
  }, mc.cores=mc.cores)

  prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)

  return(list(prob = prob, penalties = penalties, subsamples = subsamples))
}
