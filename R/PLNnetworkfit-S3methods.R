## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNnetworkfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNnetworkfit <- function(Robject) {inherits(Robject, "PLNnetworkfit")}

#' Extract and plot the network (partial correlation, support or inverse covariance) from a [`PLNnetworkfit`] object
#'
#' @name plot.PLNnetworkfit
#'
#' @param x an R6 object with class [`PLNnetworkfit`]
#' @param type character. Value of the weight of the edges in the network, either "partial_cor" (partial correlation) or "support" (binary). Default is `"partial_cor"`.
#' @param output the type of output used: either 'igraph' or 'corrplot'. Default is `'igraph'`.
#' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
#' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
#' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
#' @param layout an optional igraph layout. Only relevant for igraph output.
#' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
#' @param ... Not used (S3 compatibility).
#'
#' @return Send back an invisible object (igraph or Matrix, depending on the output chosen) and optionally displays a graph (via igraph or corrplot for large ones)
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' fits <- PLNnetwork(Abundance ~ 1, data = trichoptera)
#' myNet <- getBestModel(fits)
#' \dontrun{
#' plot(myNet)
#' }
#' @export
plot.PLNnetworkfit <-
  function(x,
           type            = c("partial_cor", "support"),
           output          = c("igraph", "corrplot"),
           edge.color      = c("#F8766D", "#00BFC4"),
           remove.isolated = FALSE,
           node.labels     = NULL,
           layout          = layout_in_circle,
           plot            = TRUE, ...) {
    invisible(x$plot_network(type, output, edge.color, remove.isolated, node.labels, layout, plot))
  }
