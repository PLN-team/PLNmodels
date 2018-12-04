#' An R6 Class to represent a PLNfit in a LDA framework
#'
#' @description The function \code{\link{PLNLDA}} produces an instance of an object with class \code{PLNPLDAfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by  \code{\link[=PLNfit]{PLNfit}}, the
#' \code{\link[=plot.PLNLDAfit]{plot.PLNPCAfit}} method for LDA vizualization and
#' \code{\link[=predict.PLNLDAfit]{predict.PLNPCAfit}} method for prediction
#'
#' @field rank the dimension of the current model
#' @field model_par a list with the matrices associated with the estimated parameters of the PLN model: Theta (covariates), Sigma (latent covariance), B (latent loadings), P (latent position) and Mu (group means)
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
      initialize = function(grouping, responses, covariates, offsets, weights, model, control) {
        super$initialize(responses, covariates, offsets, weights, model, control)
        private$grouping <- grouping
        super$optimize(responses, covariates, offsets, weights, control)
      },
      optimize = function(X, covar, design_group, control) {
        ## extract group means
        if (ncol(covar) > 0) {
          proj_orth_X <- (diag(self$n) - covar %*% solve(crossprod(covar)) %*% t(covar))
          P <- proj_orth_X %*% self$latent_pos(X, matrix(0, self$n, self$q))
          Mu <- t(rowsum(P, private$grouping) / tabulate(private$grouping))
        } else {
          Mu <- private$Theta
        }
        colnames(Mu) <- colnames(design_group)
        private$Mu <- Mu
        nk <- table(private$grouping)
        Mu_bar <- as.vector(Mu %*% nk / self$n)
        private$B <- Mu %*% diag(nk) %*% t(Mu) / self$n - Mu_bar %o% Mu_bar
      },
      setVisualization = function(scale.unit = FALSE) {
        Wm1B <- solve(private$Sigma) %*% private$B
        private$svdLDA <- svd(scale(Wm1B,TRUE, scale.unit), nv = self$rank)
        ## P <- self$latent_pos(model.matrix( ~ private$grouping + 0), matrix(0, self$n, self$q))
        P <- private$M + tcrossprod(model.matrix( ~ private$grouping + 0), private$Mu) ## P = M + G Mu
        private$P <- scale(P, TRUE, FALSE)
      }
    ),
    private = list(
      B        = NULL,
      P        = NULL,
      Mu       = NULL,
      grouping = NULL,
      svdLDA   = NULL
    ),
    active = list(
      rank = function() {nlevels(private$grouping) - 1},
      degrees_freedom = function() {self$p * (self$d + self$rank)},
      model_par = function() {
        par <- super$model_par
        par$B  <- private$B
        par$P  <- private$P
        par$Mu <- private$Mu
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

PLNLDAfit$set("public", "postTreatment",
function(responses, covariates, offsets) {
  super$postTreatment(responses, covariates, offsets)
  rownames(private$B) <- colnames(private$B) <- colnames(responses)
  if (private$covariance != "spherical") colnames(private$S) <- 1:self$q
  self$setVisualization()
})

# Plot the individual map of a specified axis for a \code{PLNLDAfit} object
PLNLDAfit$set("public", "plot_individual_map",
  function(axes = 1:min(2,self$rank), main = "Individual Factor Map", plot = TRUE) {

    .scores <- data.frame(self$scores[,axes, drop = FALSE])
    colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
    .scores$labels <- private$grouping
    .scores$names <- rownames(private$M)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_ind_map(.scores, axes_label, main)
    if (plot) print(p)
    invisible(p)
})

# Plot the correlation circle of a specified axis for a \code{PLNLDAfit} object
PLNLDAfit$set("public", "plot_correlation_circle",
  function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "default", plot=TRUE) {

    ## data frame with correlations between variables and PCs
    correlations <- as.data.frame(self$corr_circle[, axes, drop = FALSE])
    colnames(correlations) <- paste0("axe", 1:length(axes))
    correlations$labels <- cols
    correlations$names  <- rownames(correlations)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_corr_circle(correlations, axes_label, main)

    if (plot) print(p)
    invisible(p)
})

# Plot a summary of the current \code{PLNLDAfit} object
PLNLDAfit$set("public", "plot_LDA",
  function(nb_axes = min(3, self$rank), var_cols = "default", plot = TRUE) {

    axes <- 1:nb_axes
    if (nb_axes > 1) {
      pairs.axes <- combn(axes, 2, simplify = FALSE)

      ## get back all individual maps
      ind.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_individual_map(axes = pair, plot = FALSE, main="") + theme(legend.position="none")
        return(ggplotGrob(ggobj))
      })

      ## get back all correlation circle
      cor.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_correlation_circle(axes = pair, plot = FALSE, main = "", cols = var_cols)
        return(ggplotGrob(ggobj))
      })

      ## plot that appear on the diagonal
      crit <- setNames(c(NA,NA,NA), c("loglikelihood", "BIC", "ICL"))
      criteria.text <- paste("Model Selection\n\n", paste(names(crit), round(crit, 2), sep=" = ", collapse="\n"))
      percentV.text <- paste("Axes contribution\n\n", paste(paste("axis",axes), paste0(": ", round(100*self$percent_var[axes],3), "%"), collapse="\n"))

      diag.grobs <- list(textGrob(percentV.text),
                         g_legend(self$plot_individual_map(plot=FALSE) + guides(colour = guide_legend(nrow = 4, title="classification"))),
                         textGrob(criteria.text))
      if (nb_axes > 3)
        diag.grobs <- c(diag.grobs, rep(list(nullGrob()), nb_axes - 3))


      grobs <- vector("list", nb_axes^2)
      i.cor <- 1; i.ind <- 1; i.dia <- 1
      ind <- 0
      for (i in 1:nb_axes) {
        for (j in 1:nb_axes) {
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
      p <- arrangeGrob(grobs = grobs, ncol = nb_axes)
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

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

PLNLDAfit$set("public", "predict",
              function(newdata, newOffsets, newCounts, type = c("posterior", "response"), control = list()) {
                type = match.arg(type)
                ## Problem dimensions
                n.new <- nrow(newCounts); groups <- levels(private$grouping); K <- length(groups)
### TODO: use model formula as in PLNfit.predict
                ## Compute conditional log-likelihoods of new data, using previously estimated parameters
                cond.log.lik <- matrix(0, n.new, K)
                for (k in 1:K) { ## One VE-step to estimate the conditional (variational) likelihood of each group
                  grouping <- factor(rep(groups[k], n.new), levels = groups)
                  X <- cbind(newdata, model.matrix( ~ grouping + 0))
                  cond.log.lik[, k] <- self$VEstep(X, newOffsets, newCounts, control = control)$log.lik
                }
                ## Compute posterior probabilities
                log.prior <- rep(1, n.new) %o% log( table(private$grouping) / self$n)
                log.posterior <- cond.log.lik + log.prior
                ## format output
                res <- log.posterior
                rownames(res) <- rownames(newdata)
                colnames(res) <- groups
                if (type == "response") {
                  res <- apply(res, 1, which.max)
                  res <- factor(groups[res], levels = groups)
                  names(res) <- rownames(newdata)
                }
                return(res)
              }
)

PLNLDAfit$set("public", "show",
function() {
  super$show(paste0("Linear Discriminant Analysis for Poisson Lognormal distribution\n"))
  cat("* Additional fields for LDA\n")
  cat("    $percent_var, $corr_circle, $scores \n")
  cat("* Additional S3 methods for LDA\n")
  cat("    plot.PLNLDAfit(), predict.PLNLDAfit()\n")
})

