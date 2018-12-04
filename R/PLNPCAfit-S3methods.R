## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNPCAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNPCAfit <- function(Robject) {inherits(Robject, "PLNPCAfit"       )}

#' PCA vizualiation (individual and/or variable factor map(s)) for a \code{PLNPCAfit} object
#'
#' @name plot.PLNPCAfit
#'
#' @param x an R6 object with class PLNPCAfit
#' @param map the type of output for the PCA vizualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param ind_cols a character, factor or numeric to define the color associated with the individuals. Default is "gray65"
#' @param var_cols a character, factor or numeric to define the color associated with the variables. Default is "gray65"
#' @param plot logical. Should the plot be displayed or sent back as ggplot object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
#' @param ... Not used (S3 compatibility).
#'
#' @return displays an individual and/or variable factor maps for the corresponding axes, and/or sends back a ggplot2 or gtable object
#'
#' @export
plot.PLNPCAfit <-
  function(x,
           map      = c("both", "individual", "variable"),
           nb_axes  = min(3, x$rank),
           axes     = seq.int(min(2,x$rank)),
           ind_cols = "gray65",
           var_cols = "gray65",
           plot     = TRUE,
           main = NULL,
           ...) {

    stopifnot(isPLNPCAfit(x))

    map <- match.arg(map)

    if (is.null(main)) {
      if (map == "individual")
        main <- "Individual factor map"
      if (map == "variable")
        main <- "Variable factor map"
    }

    if (map == "individual")
      p <- x$plot_individual_map(cols = ind_cols, main = main, plot = plot)
    if (map == "variable")
      p <- x$plot_correlation_circle(cols = var_cols, main = main, plot = plot)
    if (map == "both")
      p <- x$plot_PCA(nb_axes = nb_axes, ind_cols = ind_cols, var_cols = var_cols, plot = plot)

    invisible(p)
}
