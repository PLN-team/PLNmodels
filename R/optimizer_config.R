available_algorithms_nlopt <- c("MMA", "CCSAQ", "LBFGS", "LBFGS_NOCEDAL", "VAR1", "VAR2")

## -----------------------------------------------------------------
##  Series of setter to default parameters for user's main functions
##
## should be ready to pass to either nlopt or torch optimizer
PLN_param <- function(control, n, p) {
  covariance  <- ifelse(is.null(control$covariance) , "full"    , control$covariance)
  covariance  <- ifelse(is.null(control$inception)  , covariance, control$inception$vcov_model)
  backend     <- ifelse(is.null(control$backend), "nlopt", control$backend)
  ctrl <- list(
    "backend"     = backend  ,
    "algorithm"   = "CCSAQ"  , # relevant for nlopt
    "maxeval"     = 10000    , # common to nlopt and torch
    "ftol_rel"    = 1e-8     , # common to nlopt and torch
    "xtol_rel"    = 1e-6     , # common to nlopt and torch
    "maxtime"     = -1       , # relevant for nlopt (disabled)
    "ftol_abs"    = 0.0      , # relevant for nlopt (disabled)
    "xtol_abs"    = 0.0      , # relevant for nlopt (disabled)
    "learning_rate" = 0.1    , # relevant for torch
    "trace"       = 1        ,
    "covariance"  = covariance,
    # "corr_matrix" = diag(x = 1, nrow = p, ncol = p),
    # "prec_matrix" = diag(x = 1, nrow = p, ncol = p),
    "inception"   = NULL
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

PLNPCA_param <- function(control) {
  ctrl <- list(
    "algorithm"   = "CCSAQ" ,
    "maxtime"     = -1      ,
    "xtol_abs"    = -1       ,
    "ftol_abs"    = -1       ,
    "ftol_rel"    = 1e-8    ,
    "xtol_rel"    = 1e-6    ,
    "maxeval"     = 10000   ,
    "trace"       = 1       ,
    "covariance"  = "rank"
  )
  ctrl[names(control)] <- control
  stopifnot(ctrl$algorithm %in% available_algorithms_nlopt)
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
