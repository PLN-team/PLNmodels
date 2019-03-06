#' An R6 Class to represent a PLNfit in a PCA framework
#'
#' @description The function \code{\link{PLNPCA}} produces a collection of models which are instances of object with class \code{PLNPCAfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by  \code{\link[=PLNfit]{PLNfit}} and the
#' \code{\link[=plot.PLNPCAfit]{plot.PLNPCAfit}} methods for PCA vizualization
#'
#' @field rank the dimension of the current model
#' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: Theta (covariates), Sigma (latent covariance) and B (latent loadings)
#' @field var_par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field latent a matrix: values of the latent vector (Z in the model)
#' @field optim_par a list with parameters useful for monitoring the optimization
#' @field loglik variational lower bound of the loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
#' @field nb_param number of parameters in the current PLN model
#' @field percent_var the percent of variance explained by each axis
#' @field corr_circle a matrix of correlations to plot the correlation circles
#' @field scores a matrix of scores to plot the individual factor maps (a.k.a. principal comonents)
#' @field rotation a matrix of rotation of the latent space
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#' myPCA <- getBestModel(myPCAs)
#' class(myPCA)
#' print(myPCA)
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfamily]{PLNPCAfamily}}
PLNPCAfit <-
  R6Class(classname = "PLNPCAfit",
    inherit = PLNfit,
    public  = list(
      initialize = function(rank, responses, covariates, offsets, weights, model, control) {
        super$initialize(responses, covariates, offsets, weights, model, control)
        if (!is.null(control$svdM)) {
          svdM <- control$svdM
        } else {
          svdM <- svd(private$M, nu = rank, nv = self$p)
        }
        private$M <- svdM$u[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank) %*% t(svdM$v[1:rank, 1:rank, drop = FALSE])
        private$S <- matrix(max(control$lower_bound), self$n, rank)
        private$B <- svdM$v[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank)/sqrt(self$n)
        private$covariance <- "rank"
      },
      update = function(Theta=NA, Sigma=NA, B=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
        super$update(Theta = Theta, Sigma = Sigma, M = M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
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
      nb_param = function() {self$p * (self$d + self$rank)},
      entropy  = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(private$S)))},
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
        corr <- t(t(private$svdBM$v[, 1:self$rank, drop = FALSE]) * private$svdBM$d[1:self$rank]^2)
        corr <- corr/sqrt(rowSums(corr^2))
        rownames(corr) <- rownames(private$Sigma)
        corr
      },
      scores     = function() {
        scores <- t(t(private$svdBM$u[, 1:self$rank]) * private$svdBM$d[1:self$rank])
        rownames(scores) <- rownames(private$M)
        scores
      },
      rotation   = function() {
        rotation <- private$svdBM$v[, 1:self$rank, drop = FALSE]
        rownames(rotation) <- rownames(private$Sigma)
        rotation
      }
### Why not sending back the rotation matrix ?
    )
)

## Call to the C++ optimizer and update of the relevant fields
PLNPCAfit$set("public", "optimize",
function(responses, covariates, offsets, weights, control) {

  ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
  opts <- control
  opts$xtol_abs <- c(rep(0, self$p*(self$d + self$q) + self$n * self$q),
                     rep(control$xtol_abs, self$n*self$q))
  opts$lower_bound <- c(rep(-Inf, self$p*self$d), # Theta
                        rep(-Inf, self$p*self$q) , # B
                        rep(-Inf, self$n*self$q) , # M
                        rep(control$lower_bound,self$n*self$q)) # S
  opts$rank <- self$q
  optim_out <- optim_rank(
    c(private$Theta, private$B, private$M, private$S),
    responses, covariates, offsets, weights, opts
  )

  self$update(
    B     = optim_out$B,
    Theta = optim_out$Theta,
    Sigma = optim_out$Sigma,
    M     = optim_out$M,
    S     = optim_out$S,
    A     = optim_out$A,
    Z     = optim_out$Z,
    Ji    = optim_out$loglik,
    monitoring = list(
      iterations = optim_out$iterations,
      status     = optim_out$status,
      message    = statusToMessage(optim_out$status))
  )
})

PLNPCAfit$set("public", "postTreatment",
function(responses, covariates, offsets, weights, nullModel) {
  super$postTreatment(responses, covariates, offsets, weights, nullModel = nullModel)
  colnames(private$B) <- colnames(private$M) <- 1:self$q
  rownames(private$B) <- colnames(responses)
  if (private$covariance != "spherical") colnames(private$S) <- 1:self$q
  self$setVisualization()
})

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

PLNPCAfit$set("public", "plot_individual_map",
  function(axes=1:min(2,self$rank), main = "Individual Factor Map", plot = TRUE, cols = "default") {

    .scores <- data.frame(self$scores[,axes, drop = FALSE])
    colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
    .scores$labels <- cols
    .scores$names <- rownames(private$M)
    axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

    p <- get_ggplot_ind_map(.scores, axes_label, main)
    if (plot) print(p)
    invisible(p)
})

PLNPCAfit$set("public", "plot_correlation_circle",
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

#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob nullGrob
PLNPCAfit$set("public", "plot_PCA",
  function(nb_axes = min(3, self$rank), ind_cols = "ind_cols", var_cols = "var_cols", plot = TRUE) {

    axes <- 1:nb_axes
    if (nb_axes > 1) {
      pairs.axes <- combn(axes, 2, simplify = FALSE)

      ## get back all individual maps
      ind.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_individual_map(axes = pair, plot = FALSE, main = "", cols = ind_cols) + theme(legend.position = "none")
        return(ggplotGrob(ggobj))
      })

      ## get back all correlation circle
      cor.plot <- lapply(pairs.axes, function(pair) {
        ggobj <- self$plot_correlation_circle(axes = pair, plot = FALSE, main = "", cols = var_cols)
        return(ggplotGrob(ggobj))
      })

      ## plot that appear on the diagonal
      crit <- setNames(c(self$loglik, self$BIC, self$ICL), c("loglik", "BIC", "ICL"))
      criteria.text <- paste("Model Selection\n\n", paste(names(crit), round(crit, 2), sep=" = ", collapse="\n"))
      percentV.text <- paste("Axes contribution\n\n", paste(paste("axis",axes), paste0(": ", round(100*self$percent_var[axes],3), "%"), collapse="\n"))

      diag.grobs <- list(textGrob(percentV.text),
                         g_legend(self$plot_individual_map(plot=FALSE, cols=ind_cols) + guides(colour = guide_legend(nrow = 4, title="classification"))),
                         textGrob(criteria.text))
      if (nb_axes > 3)
        diag.grobs <- c(diag.grobs, rep(list(nullGrob()), nb_axes - 3))


      grobs <- vector("list", nb_axes^2)
      i.cor <- 1; i.ind <- 1; i.dia <- 1
      ind <- 0
      for (i in 1:nb_axes) {
        for (j in 1:nb_axes) {
          ind <- ind + 1
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

# Compute the (one-data) Fisher information matrix of Theta using one of two approximations scheme.
PLNPCAfit$set("public", "compute_fisher",
  function(type = c("wald", "louis"), X = NULL) {
    type = match.arg(type)
    if (type == "louis") {
      stop("Louis approximation scheme not available yet for object of class PLNPLCA, use type = \"wald\" instead.")
    }
    super$compute_fisher(type = "wald", X = X)
  }
)

PLNPCAfit$set("public", "show",
function() {
  super$show(paste0("Poisson Lognormal with rank constrained for PCA (rank = ",self$rank,")\n"))
  cat("* Additional fields for PCA\n")
  cat("    $percent_var, $corr_circle, $scores, $rotation\n")
  cat("* Additional S3 methods for PCA\n")
  cat("    plot.PLNPCAfit()\n")
})
