available_algorithms_nlopt <- c("MMA", "CCSAQ", "LBFGS", "LBFGS_NOCEDAL", "VAR1", "VAR2")

config_default_nlopt <-
  list(
    algorithm     = "CCSAQ"  ,
    maxeval       = 10000    ,
    ftol_rel      = 1e-8     ,
    xtol_rel      = 1e-6     ,
    ftol_abs      = 0.0      ,
    xtol_abs      = 0.0      ,
    maxtime       = -1
  )

config_default_torch <-
  list(
    maxeval       = 10000    ,
    ftol_rel      = 1e-8     ,
    xtol_rel      = 1e-6     ,
    learning_rate = 0.1      ,
    trace         = 1
  )

## -----------------------------------------------------------------
##  Series of setter to default parameters for user's main functions
##
## should be ready to pass to either nlopt or torch optimizer

#' Control of PLN fit
#'
#' Helper to define list of parameters to control the PLN fit. All arguments have defaults.
#'
#' @param backend optimization back used, either "nlopt" or "torch". Default is "nlopt"
#' @param covariance character setting the model for the covariance matrix. Either "full", "diagonal", "spherical", "fixed" or "genetic". Default is "full".
#' @param Omega precision matrix of the latent variables. Inverse of Sigma. Must be specified if `covariance` is "fixed"
#' @param config_optim a list for controlling the optimizer (either "nlopt" or "torch" backend). See details
#' @param trace a integer for verbosity.
#' @param inception Set up the parameters initialization: by default, the model is initialized with a multivariate linear model applied on
#'    log-transformed data, and with the same formula as the one provided by the user. However, the user can provide a PLNfit (typically obtained from a previous fit),
#'    which sometimes speeds up the inference.
#'
#' @return list of parameters configuring the fit.
#'
#' @details The list of parameters `config_optim` controls the optimizers. When "nlopt" is chosen the following entries are relevant
#' * "algorithm" the optimization method used by NLOPT among LD type, e.g. "CCSAQ", "MMA", "LBFGS". See NLOPT documentation for further details. Default is "CCSAQ".
#' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
#' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-8
#' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-6
#' * "ftol_abs" stop when an optimization step changes the objective function by less than ftol_abs. Default is 0.0 (disabled)
#' * "xtol_abs" stop when an optimization step changes every parameters by less than xtol_abs. Default is 0.0 (disabled)
#' * "maxtime" stop when the optimization time (in seconds) exceeds maxtime. Default is -1 (disabled)
#'
#' When "torch" backend is used, with the following entries are relevant:
#' * "maxeval" stop when the number of iteration exceeds maxeval. Default is 10000
#' * "ftol_rel" stop when an optimization step changes the objective function by less than ftol multiplied by the absolute value of the parameter. Default is 1e-8
#' * "xtol_rel" stop when an optimization step changes every parameters by less than xtol multiplied by the absolute value of the parameter. Default is 1e-6
#'
#' @export
PLN_param <- function(
    backend       = "nlopt",
    trace         = 1      ,
    covariance    = "full" ,
    Omega         = NULL   ,
    config_optim  = list() ,
    inception     = NULL     # pretrained PLNfit used as initialization
) {
  stopifnot(backend %in% c("nlopt", "torch"))
  stopifnot(config_optim$algorithm %in% available_algorithms_nlopt)
  if (covariance == "fixed") stopifnot(inherits(Omega, "matrix"))
  if (backend == "nlopt") {
    config <- config_default_nlopt
    config[names(config_optim)] <- config_optim
  }
  if (backend == "torch") {
    config <- config_default_torch
    config[names(config_optim)] <- config_optim
  }
  if(!is.null(inception)) stopifnot(isPLNfit(inception))
  structure(list(
    backend       = backend   ,
    trace         = trace     ,
    covariance    = covariance,
    Omega         = Omega     ,
    config_optim  = config    ,
    covariance    = covariance,
    inception     = inception   ), class = "PLN_param")
}

PLNPCA_param <- function(control) {
  ctrl <- list(
    "backend"       = "nlopt"  ,
    "algorithm"     = "CCSAQ"  , # relevant for nlopt
    "maxeval"       = 10000    , # common to nlopt and torch
    "ftol_rel"      = 1e-8     , # common to nlopt and torch
    "xtol_rel"      = 1e-6     , # common to nlopt and torch
    "maxtime"       = -1       , # relevant for nlopt (disabled)
    "ftol_abs"      = 0.0      , # relevant for nlopt (disabled)
    "xtol_abs"      = 0.0      , # relevant for nlopt (disabled)
    "learning_rate" = 0.1      , # relevant for torch
    "trace"         = 1        ,
    "covariance"    = "rank"
  )
  ctrl[names(control)] <- control
  stopifnot(ctrl$algorithm %in% available_algorithms_nlopt)
  ctrl
}

## should be ready to pass to nlopt optimizer
PLNmixture_param <- function(control, n, p) {
  covariance  <- ifelse(is.null(control$covariance) , "spherical", control$covariance)
  covariance  <- ifelse(is.null(control$inception), covariance  , control$inception$model)
  ctrl <- list(
    "ftol_out"    = 1e-3,
    "maxit_out"   = 50,
    "algorithm"   = "CCSAQ",
    "maxeval"     = 10000  ,
    "maxtime"     = -1     ,
    "ftol_rel"    = ifelse(n < 1.5*p, 1e-6, 1e-8),
    "ftol_abs"    = 0,
    "xtol_rel"    = 1e-6,
    "xtol_abs"    = 0.0,
    "backend"     = "nlopt",
    "learning_rate" = 0.1,
    "trace"       = 1,
    "covariance"  = covariance,
    "iterates"    = 2,
    "smoothing"   = 'both',
    "inception"   = NULL,
    "init_cl"     = 'kmeans'
  )
  ctrl[names(control)] <- control
  ctrl
}

PLNnetwork_param <- function(control, n, p) {
  xtol_abs    <- ifelse(is.null(control$xtol_abs)   , 0, control$xtol_abs)
  ctrl <-  list(
    "ftol_out"  = 1e-5,
    "maxit_out" = 20,
    "warm"        = FALSE,
    "algorithm"   = "CCSAQ",
    "ftol_rel"    = ifelse(n < 1.5*p, 1e-6, 1e-8),
    "ftol_abs"    = 0       ,
    "xtol_rel"    = 1e-6    ,
    "xtol_abs"    = xtol_abs,
    "maxeval"     = 10000   ,
    "maxtime"     = -1      ,
    "trace"       = 1       ,
    "covariance"  = "sparse"
  )
  ctrl[names(control)] <- control
  stopifnot(ctrl$algorithm %in% available_algorithms_nlopt)
  ctrl
}

status_to_message <- function(status) {
  message <- switch(as.character(status),
        "1"  = "success",
        "2"  = "success, stopval was reached",
        "3"  = "success, ftol_rel or ftol_abs was reached",
        "4"  = "success, xtol_rel or xtol_abs was reached",
        "5"  = "success, maxeval was reached",
        "6"  = "success, maxtime was reached",
        "-1" = "failure",
        "-2" = "invalid arguments",
        "-3" = "out of memory.",
        "-4" = "roundoff errors led to a breakdown of the optimization algorithm",
        "-5" = "forced termination:",
        "Return status not recognized"
  )
  message
}

status_to_message_nlopt <- function(status) {
  message <- switch(as.character(status),
                    "1"  = "success",
                    "2"  = "success, stopval was reached",
                    "3"  = "success, ftol_rel or ftol_abs was reached",
                    "4"  = "success, xtol_rel or xtol_abs was reached",
                    "5"  = "success, maxeval was reached",
                    "6"  = "success, maxtime was reached",
                    "-1" = "failure",
                    "-2" = "invalid arguments",
                    "-3" = "out of memory.",
                    "-4" = "roundoff errors led to a breakdown of the optimization algorithm",
                    "-5" = "forced termination:",
                    "Return status not recognized"
  )
  message
}
