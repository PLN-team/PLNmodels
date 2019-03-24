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
#' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
#' @field nb_param number of parameters in the current PLN model
#' @field percent_var the percent of variance explained by each axis
#' @field corr_map a matrix of correlations to plot the correlation circles
#' @field scores a matrix of scores to plot the individual factor maps
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
#' class(myPLNLDA)
#' print(myPLNLDA)
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
        rownames(Mu) <- rownames(private$Theta)
        private$Mu <- Mu
        nk <- table(private$grouping)
        Mu_bar <- as.vector(Mu %*% nk / self$n)
        private$B <- Mu %*% diag(nk) %*% t(Mu) / self$n - Mu_bar %o% Mu_bar
      },
      setVisualization = function(scale.unit = FALSE) {
        Wm1B <- solve(private$Sigma) %*% private$B
        private$svdLDA <- svd(scale(Wm1B,TRUE, scale.unit), nv = self$rank)
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
      nb_param = function() {self$p * (self$d + self$rank)},
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
      corr_map = function() {
        corr <- cor(private$P, self$scores)
        rownames(corr) <- rownames(private$B)
        corr
      },
      scores     = function() {
        scores <- private$P %*% t(t(private$svdLDA$u[, 1:self$rank]) * private$svdLDA$d[1:self$rank])
        rownames(scores) <- rownames(private$M)
        scores
      },
      group_means = function() {
        self$model_par$Mu
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
PLNLDAfit$set("public", "plot_correlation_map",
  function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "default", plot=TRUE) {

    ## data frame with correlations between variables and PCs
    correlations <- as.data.frame(self$corr_map[, axes, drop = FALSE])
    colnames(correlations) <- paste0("axe", 1:length(axes))
    correlations$labels <- cols
    correlations$names  <- rownames(correlations)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_corr_square(correlations, axes_label, main)

    if (plot) print(p)
    invisible(p)
})


# Plot a summary of the current \code{PLNLDAfit} object
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid nullGrob textGrob
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
        ggobj <- self$plot_correlation_map(axes = pair, plot = FALSE, main = "", cols = var_cols)
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
        self$plot_correlation_map(plot = FALSE)
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
  function(newdata,
           type = c("posterior", "response", "scores"),
           scale = c("log", "prob"),
           prior = NULL,
           control = list(), envir = parent.frame()) {

    type = match.arg(type)

    if (type == "scores") scale <- "prob"
    scale = match.arg(scale)

    args <- extract_model(call("PLNLDA", formula = private$model, data = newdata), envir)

    ## Problem dimensions
    n.new  <- nrow(args$Y)
    p      <- ncol(args$Y)
    groups <- levels(private$grouping)
    K <- length(groups)

    ## Initialize priors
    if (is.null(prior)) {
      prior <- table(private$grouping)
    } else {
      names(prior) <- groups
    }
    if (any(prior <= 0) || anyNA(prior)) stop("Prior group proportions should be positive.")
    prior <- prior / sum(prior)

    ## Compute conditional log-likelihoods of new data, using previously estimated parameters
    cond.log.lik <- matrix(0, n.new, K)
    if (type == "scores") latent_pos <- array(0, dim = c(n.new, K, p))
    for (k in 1:K) { ## One VE-step to estimate the conditional (variational) likelihood of each group
      grouping <- factor(rep(groups[k], n.new), levels = groups)
      X <- cbind(args$X, model.matrix( ~ grouping + 0))
      # - remove intercept so that design matrix is compatible with the one used for inference
      xint <- match("(Intercept)", colnames(X), nomatch = 0L)
      if (xint > 0L) X <- X[, -xint, drop = FALSE]
      ve_step <- self$VEstep(X, args$O, args$Y, control = control)
      cond.log.lik[, k] <- ve_step$log.lik
      if (type == "scores") {
        latent_pos[ , k, ] <- ve_step$M + rep(1, n.new) %o% self$group_means[, k]
      }
    }
    ## Compute (unnormalized) posterior probabilities
    log.prior <- rep(1, n.new) %o% log(prior)
    log.posterior <- cond.log.lik + log.prior

    ## format output
    res <- log.posterior
    if (scale == "prob") { ## change log-likelihoods to probabilities
      ## trick to avoid rounding errors before exponentiation
      row_max <- apply(res, 1, max)
      res <- exp(sweep(res, 1, row_max, "-"))
      res <- sweep(res, 1, rowSums(res), "/")
    }
    rownames(res) <- rownames(newdata)
    colnames(res) <- groups

    if (type == "scores") { ## compute average latent positions
      weights <- res %o% rep(1, p)
      ## weighted average over all groups (second dimension of the array)
      average_positions <- apply(latent_pos * weights, c(1, 3), sum)
      ## Center positions using the same origin as for the learning set
      centers <- attr(private$P, 'scaled:center')
      centered_positions <- sweep(average_positions, 2, centers)
      ## transformation through svd to project in the same space as learning samples
      scores <- centered_positions %*% t(t(private$svdLDA$u[, 1:self$rank]) * private$svdLDA$d[1:self$rank])
      rownames(scores) <- rownames(private$M)
      return(scores)
    }

    if (type == "response") {
      res <- apply(res, 1, which.max)
      res <- factor(groups[res], levels = groups)
      names(res) <- rownames(newdata)
    }
    res
  }
)

PLNLDAfit$set("public", "show",
function() {
  super$show(paste0("Linear Discriminant Analysis for Poisson Lognormal distribution\n"))
  cat("* Additional fields for LDA\n")
  cat("    $percent_var, $corr_map, $scores, $group_means\n")
  cat("* Additional S3 methods for LDA\n")
  cat("    plot.PLNLDAfit(), predict.PLNLDAfit()\n")
})

