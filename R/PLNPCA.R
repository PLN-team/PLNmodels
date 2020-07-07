##' Poisson lognormal model towards Principal Component Analysis
##'
##' Fit the PCA variants of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
##'
##' @param formula an object of class "formula": a symbolic description of the model to be fitted.
##' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
##' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
##' @param weights an optional vector of observation weights to be used in the fitting process.
##' @param ranks a vector of integer containing the successive ranks (or number of axes to be considered)
##' @param control_init a list for controlling the optimization at initialization.  See details of function [PLN()].
##' @param control_main a list for controlling the main optimization process. See details.
##'
##' @return an R6 object with class [`PLNPCAfamily`], which contains
##' a collection of models with class [`PLNPCAfit`]
##'
##' @details The list of parameters `control_main` controls the optimization of the main process, with the following entries:
##' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-8
##' * "ftol_abs" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0
##' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4
##' * "xtol_abs" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 0
##' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
##' * "maxtime" stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)
##' * "algorithm" the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS",
##'     "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".
##' * "trace" integer for verbosity. Useless when `cores` > 1
##' * "cores" The number of core used to parallelize jobs over the `ranks` vector. Default is 1.
##'
##'
##' @rdname PLNPCA
##' @examples
##' data(trichoptera)
##' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
##' myPCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
##' @seealso The classes [`PLNPCAfamily`] and [`PLNPCAfit`]
##' @importFrom stats model.frame model.matrix model.response model.offset
##' @export
PLNPCA <- function(formula, data, subset, weights, ranks = 1:5, control_init = list(), control_main = list()) {

  ## extract the data matrices and weights
  args <- extract_model(match.call(expand.dots = FALSE), parent.frame())

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl_init <- PLN_param(control_init, nrow(args$Y), ncol(args$Y), ncol(args$X))
  ctrl_main <- PLNPCA_param(control_main)

  ## Instantiate the collection of PLN models, initialized by PLN with full rank
  if (ctrl_main$trace > 0) cat("\n Initialization...")
  myPCA <- PLNPCAfamily$new(ranks, args$Y, args$X, args$O, args$w,
                            args$model, args$xlevels, ctrl_init)

  ## Adjust the PLN models
  if (ctrl_main$trace > 0) cat("\n\n Adjusting", length(ranks), "PLN models for PCA analysis.\n")
  myPCA$optimize(ctrl_main)

  ## Post-treatments: pseudo-R2, rearrange criteria and prepare PCA visualization
  if (ctrl_main$trace > 0) cat("\n Post-treatments")
  myPCA$postTreatment()

  if (ctrl_main$trace > 0) cat("\n DONE!\n")
  myPCA
}
