##' Poisson lognormal mixture model
##'
##' Fit the mixture variants of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param clusters a vector of integer containing the successive number of clusters (or components) to be considered
##' @param control_init a list for controling the optimization at initialization. See details.
##' @param control_main a list for controling the main optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNMMfamily]{PLNMMfamily}}, which contains
##' a collection of models with class \code{\link[=PLNMMfit]{PLNMMfit}}
##'
##' @details The list of parameters \code{control_init} and \code{control_main} control the optimization of the intialization and the main process, with the following entries
##' \itemize{
##'   \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol_rel multiplied by the absolute value of the parameter.}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol_abs .}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol_rel multiplied by the absolute value of the parameter.}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol_abs.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"method"}{the optimization method used by NLOPT among LD type, i.e. "CCSAQ"  "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "TNEWTON_VAR1", "TNEWTON_VAR2". See NLOPTR documentation for further details. Default is "CCSAQ".}
##'  \item{"lbvar"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-5.}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##'  \item{"inception"}{a PLNfit to start with. If NULL, a PLN is fitted on the . If an R6 object with class 'PLNfit' is given, it is used to initialize the model.}
##'  \item{"cores"}{The number of core used to paralellize jobs over the \code{ranks} vector. Default is 1.}
##' }
##'
##' @rdname PLNMM
##' @examples
##' ## See the vignette: vignette("trichoptera", package="PLNmodels")
##' @seealso The classes \code{\link[=PLNMM]{PLNMMfamily}} and \code{\link[=PLNMMfit]{PLNMMfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNMM <- function(formula, data, subset, clusters = 1:10,  control_init = list(), control_main = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())
  # be sure that the intercept is in: to avoid conflict with the cluster means
  # - remove intercept so that 'grouping' describes group means
  xint <- match("(Intercept)", colnames(args$X), nomatch = 0L)
  if (xint == 0L) args$X <- cbind(1, args$X)

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl_init <- PLNMM_param(control_init, nrow(args$Y), ncol(args$Y), "init")
  ctrl_main <- PLNMM_param(control_main, nrow(args$Y), ncol(args$Y), "main")

  ## Instantiate the collection of PLN models
  if (ctrl_main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNMMfamily$new(
    clusters   = clusters,
    responses  = args$Y,
    covariates = args$X,
    offsets    = args$O,
    control    = ctrl_init
  )

  ## Now adjust the PLN models
  if (ctrl_main$trace > 0) cat("\n\n Adjusting", length(clusters), "PLN mixture models.\n")
  myPLN$optimize(ctrl_main)

  ## Post-treatments: Compute pseudo-R2, rearrange criteria and the visualization for PCA
  # if (ctrl_main$trace > 0) cat("\n Post-treatments")
  # myPLN$postTreatment()

  if (ctrl_main$trace > 0) cat("\n DONE!\n")
  myPLN
}
