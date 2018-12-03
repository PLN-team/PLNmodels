## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNPCAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNPCAfit <- function(Robject) {inherits(Robject, "PLNPCAfit"       )}

#' Plot the individual or variable factor map for a \code{PLNPCAfit} object
#'
#' @name plot.PLNPCAfit
#'
#' @param x an R6 object with class PLNPCAfit
#' @param axes numeric, the axes to use for the plot. Default it c(1,2)
#' @param plot logical. Should the plot be displayed or sent back (ggplot object)
#'
#' @param ... additional parameters: cols, main.
#'
#' @return displays a individual or variable map for the corresponding axes and/or sends back a ggplot2 object
#'
#' @export
plot.PLNPCAfit <-
  function(x,
           map  = c("individual", "variable"),
           axes = 1:min(2,x$rank),
           plot = TRUE, ...) {

    stopifnot(isPLNPCAfit(x))

    map <- match.arg(map)

    param <- list(...)
    if (is.null(param$cols)) cols <- "gray65" else cols <- param$cols
    if (is.null(param$main)) {
      if (map == "individual")
        main <- "Individual factor map"
      if (map == "variable")
        main <- "Variable factor map"
    } else {
      main <- param$main
    }

    if (map == "individual")
      p <- x$plot_individual_map(cols = cols, main = main, plot = FALSE)
    if (map == "variable")
      p <- x$plot_correlation_circle(cols = cols, main = main, plot = FALSE)
    if (plot)
      print(p)

    invisible(p)
}

