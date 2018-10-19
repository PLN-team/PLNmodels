##' Poisson lognormal model towards Principal Component Analysis
##'
##' Fit the PCA variants of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param ranks a vector of integer containing the successive ranks (or number of axes to be considered)
##' @param control_init a list for controling the optimization at initialization.  See details of function \code{\link[=PLN]{PLN}}.
##' @param control_main a list for controling the main optimization process. See details.
##'
##' @return an R6 object with class \code{\link[=PLNPCAfamily]{PLNPCAfamily}}, which contains
##' a collection of models with class \code{\link[=PLNPCAfit]{PLPCAfit}}
##'
##' @details The list of parameters \code{control_main} controls the optimization of the main process, with the following entries
##' \itemize{
##'  \item{"ftol_rel"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-8}
##'  \item{"ftol_abs"}{stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0}
##'  \item{"xtol_rel"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"xtol_abs"}{stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4}
##'  \item{"lower_bound"}{the lower bound (box constraint) for the variational variance parameters. Default is 1e-4.}
##'  \item{"maxeval"}{stop when the number of iteration exceeds maxeval. Default is 10000}
##'  \item{"maxtime"}{stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)}
##'  \item{"algorithm"}{the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS",
##'     "TNEWTON", "TNEWTON_RESTART", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART",
##'     "TNEWTON_VAR1", "TNEWTON_VAR2". See NLOPT documentation for further details. Default is "CCSAQ".}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##'  \item{"cores"}{The number of core used to paralellize jobs over the \code{ranks} vector. Default is 1.}
##' }
##'
##' @rdname PLNPCA
##' @examples
##' data(trichoptera)
##' TotalCount <- matrix(
##' rowSums(trichoptera$Abundance),
##'   nrow = nrow(trichoptera$Abundance),
##'   ncol = ncol(trichoptera$Abundance)
##' )
##'
##' myPCA <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:6)
##' @seealso The classes \code{\link[=PLNPCAfamily]{PLNPCAfamily}} and \code{\link[=PLNPCAfit]{PLNPCAfit}}
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNPCA <- function(formula, data, subset, ranks = 1:5,  control_init = list(), control_main = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl_init <- PLN_param(control_init, nrow(args$Y), ncol(args$Y), ncol(args$X))
  if (is.null(ctrl_init$inception)) ctrl_init$inception <- ifelse(nrow(args$Y) >= 1.5*ncol(args$Y), "PLN", "LM")
  ctrl_main <- PLNPCA_param(control_main)

  ## Instantiate the collection of PLN models, initialized by PLN with full rank
  if (ctrl_main$trace > 0) cat("\n Initialization...")
  myPLN <- PLNPCAfamily$new(ranks = ranks, responses = args$Y, covariates = args$X, offsets = args$O, control = ctrl_init)

  ## Now adjust the PLN models
  if (ctrl_main$trace > 0) cat("\n\n Adjusting", length(ranks), "PLN models for PCA analysis.\n")
  myPLN$optimize(ctrl_main)

  ## Post-treatments: Compute pseudo-R2, rearrange criteria and the visualization for PCA
  if (ctrl_main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl_main$trace > 0) cat("\n DONE!\n")
  myPLN
}
