#' An R6 Class to represent a PLNfit in a network inference framework
#'
#' @description The function \code{\link{PLNnetwork}} produces a collection of models which are instances of object with class \code{PLNnetworkfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user: plot_network + methods inhererited from PLNfit.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally.
#'
#' @field penalty the level of sparsity in the current model
#' @field model.par a list with two matrices, B and Theta, which are the estimated parameters of the pPCA model
#' @field variation.par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field criteria a named vector with the value of some criteria (variational lower bound J, BIC, ICL, R2, lmin and lmax) for the different models.
#' @field convergence quantities usefull for monitoring the optimization
#' @include PLNnetworkfit-class.R
#' @importFrom R6 R6Class
#' @import Matrix
#' @import igraph
#' @importFrom corrplot corrplot
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfamily]{PLNnetworkfamily}}
PLNnetworkfit <-
  R6Class(classname = "PLNnetworkfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(penalty=NA, Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA,
                            J=NA, BIC=NA, ICL=NA, R2=NA, status=NA, iter=NA) {
        super$initialize(Theta, Sigma, Omega, M, S, J, BIC, ICL, R2, status, iter)
        private$lambda <- penalty
      },
      update = function(penalty=NA, Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA,
                        J=NA, BIC=NA, ICL=NA, R2=NA,status=NA, iter=NA) {
        super$update(Theta, Sigma, Omega, M, S, J, BIC, ICL, R2, status, iter)
        if (!anyNA(penalty)) private$lambda <- penalty
      }
    ),
    private = list(
      lambda = NULL
    ),
    active = list(
      penalty = function() {private$lambda}
    )
)

PLNnetworkfit$set("public", "latentNetwork",
  function(weighted=FALSE) {
    if (weighted) {
      res  <- abs(private$Omega)
    } else {
      res  <- 1*(private$Omega != 0)
    }
    diag(res) <- 0
    Matrix(res)
  }
)

PLNnetworkfit$set("public", "plot_network",
  function(plot=TRUE, remove.isolated = TRUE, layout=NULL) {
    net <- self$latentNetwork(weighted=TRUE)
    G <-  graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)
    if (!is.null(colnames(net)))
      V(G)$label <- colnames(net)
    weight <- (log(1+E(G)$weight)+.4) / max(log(1+E(G)$weight)+.4)
    V.deg <- degree(G)/sum(degree(G))
    V(G)$label.cex <- V.deg / max(V.deg) + .5
    V(G)$size <- V.deg * 100
    V(G)$label.color <- rgb(0, 0, .2, .8)
    V(G)$frame.color <- NA
    E(G)$color <- rgb(.5, .5, 0, weight)
    E(G)$width <- E(G)$weight * 2

    if (remove.isolated)
      G <- delete.vertices(G, which(degree(G) == 0))

    if (plot) {
      par(mfrow=c(1,2))
      par(mar=c(0.1,0.1,0.1,0.1))
      if (ncol(net) > 100)
        colnames(net) <- rownames(net) <- rep(" ", ncol(net))
      corrplot(net, method="color", is.corr=FALSE, cl.pos="n", tl.cex=0.5)
      if (!is.null(layout))
        plot(G, layout=layout)
      else
        plot(G)
      par(mfrow=c(1,1))
      par(mar=c(5.1,4.1,4.1,2.1))
    }
    invisible(G)
})
