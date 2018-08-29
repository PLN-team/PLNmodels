#' An R6 Class to represent a PLNfit in a PCA framework
#'
#' @description The function \code{\link{PLNPCA}} produces a collection of models which are instances of object with class \code{PLNPCAfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNPCAfit_plot_PCA]{plot_PCA}}, \code{\link[=PLNPCAfit_plot_individual_map]{plot_individual_map}}
#' and \code{\link[=PLNPCAfit_plot_correlation_circle]{plot_correlation_circle}}
#'
#' @field rank the dimension of the current model
#' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: Theta (covariates), Sigma (latent covariance) and B (latent loadings)
#' @field var_par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field optim_par a list with parameters useful for monitoring the optimization
#' @field loglik variational lower bound of the loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field criteria a vector with loglik, BIC, ICL, R_squared and degrees of freedom
#' @field degrees_freedom number of parameters in the current PLN model
#' @field percent_var the percent of variance explained by each axis
#' @field corr_circle a matrix of correlations to plot the correlation circles
#' @field scores a matrix of scores to plot the individual factor maps
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfamily]{PLNPCAfamily}}
PLNPCAfit <-
  R6Class(classname = "PLNPCAfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(Theta=NA, Sigma=NA, B=NA, M=NA, S=NA, J=NA, monitoring=NA) {
        super$initialize(Theta, Sigma, M, S, J, monitoring)
        private$B <- B
      },
      update = function(Theta=NA, Sigma=NA, B=NA, M=NA, S=NA, J=NA, R2=NA, monitoring=NA) {
        super$update(Theta, Sigma, M, S, J, R2, monitoring)
        if (!anyNA(B)) private$B <- B
      },
      setVisualization = function(scale.unit=FALSE) {
        P <- t(tcrossprod(private$B, private$M))
        private$svdBM <- svd(scale(P,TRUE, scale.unit), nv = self$rank)
      }
    ),
    private = list(
      B     = NULL,
      svdBM = NULL
    ),
    active = list(
      rank = function() {ncol(private$B)},
      degrees_freedom = function() {self$p * (self$d + self$rank)},
      model_par = function() {
        par <- super$model_par
        par$B <- private$B
        par
      },
      percent_var = function() {
        eigen.val <- private$svdBM$d[1:self$rank]^2
        round(eigen.val/sum(eigen.val),4)
      },
      corr_circle = function() {
        corr <- t(t(private$svdBM$v[, 1:self$rank]) * private$svdBM$d[1:self$rank]^2)
        corr <- corr/sqrt(rowSums(corr^2))
        rownames(corr) <- rownames(private$B)
        corr
      },
      scores     = function() {
        scores <- t(t(private$svdBM$u[, 1:self$rank]) * private$svdBM$d[1:self$rank])
        rownames(scores) <- rownames(private$M)
        scores
      }
      rotation   = function() {
        rotation <- private$svdBM$v[, 1:self$rank]
        rownames(rotation) <- rownames(private$B)
        rotation
      }
    )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR INTERNAL USE -> PLNfamily
## ----------------------------------------------------------------------
## Should only be accessed BY PLNfamily but R6 friend class don't exist

# Positions in the (euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#
# @name PLNfit_latent_pos
#
# @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
# @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
PLNPCAfit$set("public", "latent_pos",
function(covariates, offsets) {
  latentPos <- tcrossprod(private$M, private$B) + tcrossprod(covariates, private$Theta) + offsets
  latentPos
})

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

#' Plot the individual map of a specified axis for a \code{PLNPCAfit} object
#'
#' @name PLNPCAfit_plot_individual_map
#' @param axes numeric, the axes to use for the plot. Default it c(1,2)
#' @param main character, the title. Default is "Individual Factor Map"
#' @param plot logical. Should the plot be displayed or sent back as a ggplot object
#' @param cols a character, factor or numeric to defined the color associated with the observations. Default is "gray"
#' @return displays a individual map for thecorresponding axes and/or sends back a ggplot2 object
NULL
PLNPCAfit$set("public", "plot_individual_map",
  function(axes=1:min(2,self$rank), main="Individual Factor Map", plot=TRUE, cols="gray65") {

    .scores <- data.frame(self$scores[,axes, drop = FALSE])
    colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
    .scores$labels <- cols
    .scores$names <- rownames(private$M)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_ind_map(.scores, axes_label, main)
    if (plot) print(p)
    invisible(p)
})

#' Plot the correlation circle of a specified axis for a \code{PLNPCAfit} object
#'
#' @name PLNPCAfit_plot_correlation_circle
#' @param axes numeric, the axes to use for the plot. Default it c(1,2)
#' @param main character, the title. Default is "Variable Factor map"
#' @param plot logical. Should the plot be displayed or sent back (ggplot object)
#' @param cols a character, factor or numeric to defined the color associated with the variable. Default is "gray"
#' @return displays a correlation circle for the corresponding axes and/or sends back a ggplot2 object
NULL
PLNPCAfit$set("public", "plot_correlation_circle",
  function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "gray65", plot=TRUE) {

    ## data frame with correlations between variables and PCs
    correlations <- as.data.frame(self$corr_circle[, axes, drop = FALSE])
    colnames(correlations) <- paste0("axe", 1:length(axes))
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_corr_circle(correlations, axes_label, main, cols)

    if (plot) print(p)
    invisible(p)
})

#' Plot a summary of the current \code{PLNPCAfit} object
#'
#' @name PLNPCAfit_plot_PCA
#' @param axes numeric a vector of axes to be considered. The default is 1:min(3,rank).
#' @param plot logical. Should the plot be displayed or sent back (ggplot object)
#' @param cols.ind a character, factor or numeric to define the color associated with the individuals. Default is "gray"
#' @param var.cols a character, factor or numeric to define the color associated with the variables. Default is "gray"
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob
#' @return a plot with a matrix-like layout with size nb.axes x nb.axes, displaying individual maps and correlation circles for the corresponding axes
NULL
PLNPCAfit$set("public", "plot_PCA",
  function(cols.ind = "gray", var.cols = "gray", plot=TRUE, axes=1:min(3,self$rank)) {

    nb.axes <- length(axes)
    if (nb.axes > 1) {
      pairs.axes <- combn(axes, 2, simplify = FALSE)

      ## get back all individual maps
      ind.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_individual_map(axes=pair, plot=FALSE, main="", cols=cols.ind) + theme(legend.position="none")
        return(ggplotGrob(ggobj))
      })

      ## get back all correlation circle
      cor.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_correlation_circle(axes=pair, plot=FALSE, main="", cols = var.cols)
        return(ggplotGrob(ggobj))
      })

      ## plot that appear on the diagonal
      crit <- setNames(c(self$loglik, self$BIC, self$ICL), c("loglikelihood", "BIC", "ICL"))
      criteria.text <- paste("Model Selection\n\n", paste(names(crit), round(crit, 2), sep=" = ", collapse="\n"))
      percentV.text <- paste("Axes contribution\n\n", paste(paste("axis",axes), paste0(": ", round(100*self$percent_var[axes],3), "%"), collapse="\n"))

      diag.grobs <- list(textGrob(percentV.text, just="left"),
                         g_legend(self$plot_individual_map(plot=FALSE, cols=cols.ind) + guides(colour = guide_legend(nrow = 4, title="classification"))),
                         textGrob(criteria.text, just="left"))
      if (nb.axes > 3)
        diag.grobs <- c(diag.grobs, rep(list(nullGrob()), nb.axes-3))


      grobs <- vector("list", nb.axes^2)
      i.cor <- 1; i.ind <- 1; i.dia <- 1
      ind <- 0
      for (i in 1:nb.axes) {
        for (j in 1:nb.axes) {
          ind <- ind+1
          if (j > i) { ## upper triangular  -> cor plot
            grobs[[ind]] <- cor.plot[[i.ind]]
            i.ind <- i.ind + 1
          } else if (i == j) { ## diagonal
            grobs[[ind]] <- diag.grobs[[i.dia]]
            i.dia <- i.dia + 1
          } else {
            grobs[[ind]] <- ind.plot[[i.cor]]
            i.cor <- i.cor + 1
          }
        }
      }
      p <- arrangeGrob(grobs=grobs, ncol=nb.axes)
    } else {
      p <- arrangeGrob(grobs = list(
        self$plot_individual_map(plot = FALSE),
        self$plot_correlation_circle(plot = FALSE)
      ), ncol = 1)
    }

    if (plot)
      grid.arrange(p)

    invisible(p)
  }
)

PLNPCAfit$set("public", "show",
function() {
  super$show(paste0("Poisson Lognormal with rank constrained for PCA (rank = ",self$rank,")\n"))
  cat("* Additional fields for PCA\n")
  cat("    $percent_var, $corr_circle, $scores, $rotation \n")
  cat("* Additional methods for PCA\n")
  cat("    $plot_PCA(), $plot_correlation_circle(), $plot_individual_map() \n")
})

