## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNLDAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNLDAfit <- function(Robject) {inherits(Robject, "PLNLDAfit"       )}

#' LDA vizualiation (individual and/or variable factor map(s)) for a \code{PLNPCAfit} object
#'
#' @name plot.PLNLDAfit
#'
#' @param x an R6 object with class PLNPCAfit
#' @param map the type of output for the PCA vizualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param var_cols a character or factor to define the color associated with the variables. By default, all variables receive the default color of the current palette.
#' @param plot logical. Should the plot be displayed or sent back as ggplot object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
#' @param ... Not used (S3 compatibility).
#'
#' @return displays an individual and/or variable factor maps for the corresponding axes, and/or sends back a ggplot2 or gtable object
#'
#' @export
plot.PLNLDAfit <-
  function(x,
           map      = c("both", "individual", "variable"),
           nb_axes  = min(3, x$rank),
           axes     = seq.int(min(2,x$rank)),
           var_cols = "var_colors",
           plot     = TRUE,
           main = NULL,
           ...) {

    stopifnot(isPLNLDAfit(x))

    map <- match.arg(map)

    if (is.null(main)) {
      if (map == "individual")
        main <- "Individual factor map"
      if (map == "variable")
        main <- "Variable factor map"
    }

    if (map == "individual")
      p <- x$plot_individual_map(main = main, plot = plot)
    if (map == "variable")
      p <- x$plot_correlation_circle(cols = var_cols, main = main, plot = plot)
    if (map == "both")
      p <- x$plot_LDA(nb_axes = nb_axes, var_cols = var_cols, plot = plot)

    invisible(p)
}

#' Predict group of new samples
#'
#' @name predict.PLNLDAfit
#'
#' @param object an R6 object with class PLNLDAfit
#' @param newdata A data frame in which to look for variables with which to predict.
#' @param newOffsets A matrix in which to look for offsets with which to predict.
#' @param newCounts A matrix in which to look for counts with to predict
#' @param type The type of prediction required. The default are posterior probabilities for each group (in log-scale),
#'             the alternative "response" is the group with maximal posterior probability.
#' @param control a list for controlling the optimization. See \code{\link[=PLN]{PLN}} for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of predicted scores for each group (if type = "score") or a vector of predicted
#'         groups (if type = "response").
#' @export
predict.PLNLDAfit <- function(object, newdata,
                              type = c("posterior", "response"), control = list(), ...) {
  stopifnot(isPLNLDAfit(object))
  object$predict(newdata, type, control)
}
