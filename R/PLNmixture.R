#' Poisson lognormal mixture model
#'
#' Fit the mixture variants of the Poisson lognormal with a variational algorithm. Use the (g)lm syntax for model specification (covariates, offsets).
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param clusters a vector of integer containing the successive number of clusters (or components) to be considered
#' @param control_init a list for controlling the optimization at initialization. See details.
#' @param control_main a list for controlling the main optimization process. See details.
#'
#' @return an R6 object with class \code{\link[=PLNmixturefamily]{PLNmixturefamily}}, which contains
#' a collection of models with class \code{\link[=PLNmixturefit]{PLNmixturefit}}
#'
#' @details The list of parameters \code{control_init} and \code{control_main} control the optimization of the initialization and the main process, with the following entries
#' * "covariance" character setting the model for the covariance matrices of the mixture components. Either "full", "diagonal" or "spherical". Default is "spherical".
#' * "trace" integer for verbosity.
#' * "inception" Set up the initialization. By default, the model is initialized with a multivariate linear model applied on log-transformed data, and with the same formula as the one provided by the user. However, the user can provide a PLNfit (typically obtained from a previous fit), which sometimes speeds up the inference.
#' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-6 when n < p, 1e-8 otherwise.
#' * "ftol_abs" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 0
#' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-4
#' * "xtol_abs" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 0
#' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
#' * "maxtime" stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (no restriction)
#' * "algorithm" the optimization method used by NLOPT among LD type, i.e. "CCSAQ", "MMA", "LBFGS", "VAR1", "VAR2". See NLOPT documentation for further details. Default is "CCSAQ".
#' * "ftol_out" outer solver stops when an optimization step changes the objective function by less than xtol multiply by the absolute value of the parameter. Default is 1e-6
#' * "maxit_out" outer solver stops when the number of iteration exceeds out.maxit. Default is 50
#' * "smoothing" The smoothing to apply. Either, 'forward', 'backward' or 'both'. Default is 'both'.
#' * "iterates" number of forward/backward iteration of smoothing. Default is 2.
#'
#' @rdname PLNmixture
#' @examples
#' ## Use future to dispatch the computations on 2 workers
#' \dontrun{
#' future::plan("multisession", workers = 2)
#' }
#'
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myMixtures <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
#'                control_main = list(smoothing = "forward", iterates = 1))
#'
#' # Shut down parallel workers
#' \dontrun{
#' future::plan("sequential")
#' }
#' @seealso The classes \code{\link[=PLNmixture]{PLNmixturefamily}} and \code{\link[=PLNmixturefit]{PLNmixturefit}}
#' @importFrom stats model.frame model.matrix model.response model.offset update.formula
#' @export
PLNmixture <- function(formula, data, subset, clusters = 1:5,  control_init = list(), control_main = list()) {

  # remove the intercept term if any (will be used to deal with group means)
  the_call <- match.call(expand.dots = FALSE)
  the_call$formula <- update.formula(formula(the_call), ~ . -1)

  ## extract the data matrices and weights
  args <- extract_model(the_call, parent.frame())

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl_init <- PLN_param(control_init, nrow(args$Y), ncol(args$Y))
  ctrl_main <- PLNmixture_param(control_main, nrow(args$Y), ncol(args$Y))
  ctrl_init$covariance <- ctrl_main$covariance
  ctrl_init$init_cl    <- ctrl_main$init_cl

  ## Instantiate the collection of PLN models
  if (ctrl_main$trace > 0) cat("\n Initialization...")
  if (ctrl_main$trace > 0) cat("\n\n Adjusting", length(clusters), "PLN mixture models.\n")
  myPLN <- PLNmixturefamily$new(clusters, args$Y, args$X, args$O, args$formula, args$xlevels, ctrl_init)

  ## Now adjust the PLN models
  myPLN$optimize(ctrl_main)

  ## Smoothing to avoid local minima
  if (ctrl_main$trace > 0) cat("\n\n Smoothing PLN mixture models.\n")
  myPLN$smooth(ctrl_main)

  ## Post-treatments: Compute pseudo-R2, rearrange criteria and the visualization for PCA
  if (ctrl_main$trace > 0) cat("\n Post-treatments")
  myPLN$postTreatment()

  if (ctrl_main$trace > 0) cat("\n DONE!\n")
  myPLN
}
