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
