#' An R6 Class to represent a PLNfit in a LDA framework
#'
#' @description The function \code{\link{PLNLDA}} produces an instance of an object with class \code{PLNPLDAfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNLDAfit_plot_LDA]{plot_LDA}}, \code{\link[=PLNLDAfit_plot_individual_map]{plot_individual_map}}
#' and \code{\link[=PLNLDAfit_plot_correlation_circle]{plot_correlation_circle}}
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
#' @seealso The function \code{\link{PLNLDA}}.
PLNLDAfit <-
  R6Class(classname = "PLNLDAfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(Theta=NA, Sigma=NA, grouping = NA, M=NA, S=NA, J=NA, monitoring=NA) {
        super$initialize(Theta, Sigma, M, S, J, monitoring)
        private$grouping <- as.factor(grouping)
        nk <- table(grouping)
        mu_bar <- as.vector(Theta %*% nk / self$n)
        private$B <- Theta %*% diag(nk) %*% t(Theta) / self$n - mu_bar %o% mu_bar
      },
      setVisualization = function(scale.unit = FALSE) {
        Wm1B <- solve(private$Sigma) %*% private$B
        private$svdLDA <- svd(scale(Wm1B,TRUE, scale.unit), nv = self$rank)
        P <- self$latent_pos(model.matrix( ~ private$grouping + 0), matrix(0, self$n, self$q))
        private$P <- scale(P, TRUE, FALSE)
      }
    ),
    private = list(
      B        = NULL,
      P        = NULL,
      grouping = NULL,
      svdLDA   = NULL
    ),
    active = list(
      rank = function() {nlevels(private$grouping) - 1},
      degrees_freedom = function() {self$p * (self$d + self$rank)},
      model_par = function() {
        par <- super$model_par
        par$B <- private$B
        par$P <- private$P
        par
      },
      percent_var = function() {
        eigen.val <- private$svdLDA$d[1:self$rank]^2
        round(eigen.val/sum(eigen.val),4)
      },
      corr_circle = function() {
        corr <- cor(private$P, self$scores)
        rownames(corr) <- rownames(private$B)
        corr
      },
      scores     = function() {
        scores <- private$P %*% t(t(private$svdLDA$u[, 1:self$rank]) * private$svdLDA$d[1:self$rank])
        rownames(scores) <- rownames(private$M)
        scores
      }
    )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

#' Plot the individual map of a specified axis for a \code{PLNLDAfit} object
#'
#' @name PLNLDAfit_plot_individual_map
#' @param axes numeric, the axes to use for the plot. Default it c(1,2)
#' @param main character, the title. Default is "Individual Factor Map"
#' @param plot logical. Should the plot be displayed or sent back as a ggplot object
#' @param cols a character, factor or numeric to defined the color associated with the observations. Default is "gray"
#' @return displays a individual map for thecorresponding axes and/or sends back a ggplot2 object
NULL
PLNLDAfit$set("public", "plot_individual_map",
  function(axes=1:min(2,self$rank), main="Individual Factor Map", plot=TRUE, cols = private$grouping) {

    .scores <- data.frame(self$scores[,axes, drop = FALSE])
    colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
    .scores$labels <- cols
    .scores$names <- rownames(private$M)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_ind_map(.scores, axes_label, main)
    if (plot) print(p)
    invisible(p)
})

#' Plot the correlation circle of a specified axis for a \code{PLNLDAfit} object
#'
#' @name PLNLDAfit_plot_correlation_circle
#' @param axes numeric, the axes to use for the plot. Default it c(1,2)
#' @param main character, the title. Default is "Variable Factor map"
#' @param plot logical. Should the plot be displayed or sent back (ggplot object)
#' @param cols a character, factor or numeric to defined the color associated with the variable. Default is "gray"
#' @return displays a correlation circle for the corresponding axes and/or sends back a ggplot2 object
NULL
PLNLDAfit$set("public", "plot_correlation_circle",
  function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "gray65", plot=TRUE) {

    ## data frame with correlations between variables and PCs
    correlations <- as.data.frame(self$corr_circle[, axes, drop = FALSE])
    colnames(correlations) <- paste0("axe", 1:length(axes))
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_corr_circle(correlations, axes_label, main, cols)

    if (plot) print(p)
    invisible(p)
})

#' Plot a summary of the current \code{PLNLDAfit} object
#'
#' @name PLNLDAfit_plot_LDA
#' @param axes numeric a vector of axes to be considered. The default is 1:min(3,rank).
#' @param plot logical. Should the plot be displayed or sent back (ggplot object)
#' @param var.cols a character, factor or numeric to define the color associated with the variables. Default is "gray"
#' @return a plot with a matrix-like layout with size nb.axes x nb.axes, displaying individual maps and correlation circles for the corresponding axes
NULL
PLNLDAfit$set("public", "plot_LDA",
  function(var.cols = "gray", plot=TRUE, axes=1:min(3, self$rank)) {

    nb.axes <- length(axes)
    pairs.axes <- combn(axes, 2, simplify = FALSE)

    ## get back all individual maps
    ind.plot <- lapply(pairs.axes, function(pair) {
      ggobj <- self$plot_individual_map(axes=pair, plot=FALSE, main="") + theme(legend.position="none")
      return(ggplotGrob(ggobj))
    })

    ## get back all correlation circle
    cor.plot <- lapply(pairs.axes, function(pair) {
      ggobj <- self$plot_correlation_circle(axes=pair, plot=FALSE, main="", cols = var.cols)
      return(ggplotGrob(ggobj))
    })

    ## plot that appear on the diagonal
    crit <- setNames(c(NA,NA,NA), c("loglikelihood", "BIC", "ICL"))
    criteria.text <- paste("Model Selection\n\n", paste(names(crit), round(crit, 2), sep=" = ", collapse="\n"))
    percentV.text <- paste("Axes contribution\n\n", paste(paste("axis",axes), paste0(": ", round(100*self$percent_var[axes],3), "%"), collapse="\n"))

    diag.grobs <- list(textGrob(percentV.text, just="left"),
                       g_legend(self$plot_individual_map(plot=FALSE) + guides(colour = guide_legend(nrow = 4, title="classification"))),
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
    if (plot)
      grid.arrange(p)

    invisible(p)
  }
)

PLNLDAfit$set("public", "show",
function() {
  super$show(paste0("Linear Discriminant Analysis for Poisson Lognormal distirbution\n"))
  cat("* Additional fields for LDA\n")
  cat("    $percent_var, $corr_circle, $scores \n")
  cat("* Additional methods for LDA\n")
  cat("    $plot_LDA(), $plot_correlation_circle(), $plot_individual_map() \n")
})
PLNLDAfit$set("public", "print", function() self$show())
