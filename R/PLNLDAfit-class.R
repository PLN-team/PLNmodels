#' An R6 Class to represent a PLNfit in a LDA framework
#'
#' @description The function [PLNLDA()] produces an instance of an object with class [`PLNLDAfit`].
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by [PLNfit()], the [plot()] method for
#' LDA visualization and [predict()] method for prediction
#'
## Inheritance seems not to work for R6 classes
# @inheritParams PLNfit
# @inheritParams PLNLDAfit.predict
# @inheritParams PLNLDAfit.plot
#'
## Parameters common to many PLNLDAfit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param grouping a factor specifying the class of each observation used for discriminant analysis.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
## Parameters common to many PLNLDAfit graphical methods
#' @param map the type of output for the PCA visualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param control list controlling the optimization and the model
#' @param ind_cols a character, factor or numeric to define the color associated with the individuals. By default, all variables receive the default color of the current palette.
#' @param var_cols a character, factor or numeric to define the color associated with the variables. By default, all variables receive the default color of the current palette.
#' @param plot logical. Should the plot be displayed or sent back as ggplot object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
#'
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
#' class(myPLNLDA)
#' print(myPLNLDA)
#' }
#' @seealso The function \code{\link{PLNLDA}}.
PLNLDAfit <- R6Class(
  classname = "PLNLDAfit",
  inherit = PLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description Initialize a [`PLNLDAfit`] object
    initialize = function(grouping, responses, covariates, offsets, weights, formula, control) {
      private$grouping <- grouping
      super$initialize(responses, cbind(covariates, model.matrix( ~ grouping + 0)), offsets, weights, formula, control)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimization ----------------------
    #' @description Compute group means and axis of the LDA (noted B in the model) in the
    #' latent space, update corresponding fields
    #' @param X Abundance matrix.
    #' @param covariates design matrix. Automatically built from the covariates and the formula from the call
    #' @param config list controlling the optimization
    optimize = function(grouping, responses, covariates, offsets, weights, config) {
      super$optimize(responses, cbind(covariates, model.matrix( ~ grouping + 0)), offsets, weights, config)
      design_group <- model.matrix( ~ grouping + 0)
      ## extract group means
      if (ncol(covariates) > 0) {
        proj_orth_X <- (diag(self$n) - covariates %*% solve(crossprod(covariates)) %*% t(covariates))
        P <- proj_orth_X %*% ((cbind(covariates, design_group) %*% private$B) + private$M)
        Mu <- t(rowsum(P, private$grouping) / tabulate(private$grouping))
      } else {
        Mu <- t(private$B)
      }
      colnames(Mu) <- colnames(design_group)
      rownames(Mu) <- colnames(private$B)
      private$Mu <- Mu
      nk <- table(private$grouping)
      Mu_bar <- as.vector(Mu %*% nk / self$n)
      private$C <- Mu %*% diag(nk) %*% t(Mu) / self$n - Mu_bar %o% Mu_bar
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Post treatment --------------------
    #' @description  Update R2, fisher and std_err fields and visualization
    #' @param config list controlling the post-treatment
    postTreatment = function(grouping, responses, covariates, offsets, config) {
      covariates <- cbind(covariates, model.matrix( ~ grouping + 0))
      super$postTreatment(responses, covariates, offsets, config = config)
      rownames(private$C) <- colnames(private$C) <- colnames(responses)
      colnames(private$S) <- 1:self$q
      if (config$trace > 1) cat("\n\tCompute LD scores for visualization...")
      self$setVisualization()
    },

    #' @description Compute LDA scores in the latent space and update corresponding fields.
    #' @param scale.unit Logical. Should LDA scores be rescaled to have unit variance
    setVisualization = function(scale.unit = FALSE) {
      Wm1C <- solve(private$Sigma) %*% private$C
      private$svdLDA <- svd(scale(Wm1C,TRUE, scale.unit), nv = self$rank)
      P <- private$M + tcrossprod(model.matrix( ~ private$grouping + 0), private$Mu) ## P = M + G Mu
      private$P <- scale(P, TRUE, scale.unit)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods -----------------
    #' @description Plot the factorial map of the LDA
    # @inheritParams plot.PLNLDAfit
    #' @return a [`ggplot`] graphic
    plot_individual_map = function(axes = 1:min(2,self$rank), main = "Individual Factor Map", plot = TRUE) {

      .scores <- data.frame(self$scores[,axes, drop = FALSE])
      colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
      .scores$labels <- private$grouping
      .scores$names <- rownames(private$M)
      axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

      p <- get_ggplot_ind_map(.scores, axes_label, main)
      if (plot) print(p)
      invisible(p)
    },

    #' @description Plot the correlation circle of a specified axis for a [`PLNLDAfit`] object
    # @inheritParams plot.PLNLDAfit
    #' @param cols a character, factor or numeric to define the color associated with the variables. By default, all variables receive the default color of the current palette.
    #' @return a [`ggplot`] graphic
    plot_correlation_map = function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "default", plot=TRUE) {

      ## data frame with correlations between variables and PCs
      correlations <- as.data.frame(self$corr_map[, axes, drop = FALSE])
      colnames(correlations) <- paste0("axe", 1:length(axes))
      correlations$labels <- cols
      correlations$names  <- rownames(correlations)
      axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

      p <- get_ggplot_corr_square(correlations, axes_label, main)

      if (plot) print(p)
      invisible(p)
    },

    #' @description Plot a summary of the [`PLNLDAfit`] object
    # @inheritParams plot.PLNLDAfit
    #' @importFrom gridExtra grid.arrange arrangeGrob
    #' @importFrom grid nullGrob textGrob
    #' @return a [`grob`] object
    plot_LDA = function(nb_axes = min(3, self$rank), var_cols = "default", plot = TRUE) {

      axes <- 1:nb_axes
      if (nb_axes > 1) {
        pairs.axes <- combn(axes, 2, simplify = FALSE)

        ## get back all individual maps
        ind.plot <- lapply(pairs.axes, function(pair) {
          ggobj <- self$plot_individual_map(axes = pair, plot = FALSE, main="") + ggplot2::theme(legend.position="none")
          ggplot2::ggplotGrob(ggobj)
        })

        ## get back all correlation circle
        cor.plot <- lapply(pairs.axes, function(pair) {
          ggobj <- self$plot_correlation_map(axes = pair, plot = FALSE, main = "", cols = var_cols)
          ggplot2::ggplotGrob(ggobj)
        })

        ## plot that appear on the diagonal
        crit <- setNames(c(NA,NA,NA), c("loglikelihood", "BIC", "ICL"))
        criteria.text <- paste("Model Selection\n\n", paste(names(crit), round(crit, 2), sep=" = ", collapse="\n"))
        percentV.text <- paste("Axes contribution\n\n", paste(paste("axis",axes), paste0(": ", round(100*self$percent_var[axes],3), "%"), collapse="\n"))

        diag.grobs <- list(textGrob(percentV.text),
                           g_legend(self$plot_individual_map(plot=FALSE) + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 4, title="classification"))),
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
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Prediction methods --------------------
    #' @description Predict group of new samples
    # @inheritParams predict.PLNLDAfit
    #' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
    #' @param type The type of prediction required. The default are posterior probabilities for each group (in either unnormalized log-scale or natural probabilities, see "scale" for details), "response" is the group with maximal posterior probability and "scores" is the average score along each separation axis in the latent space, with weights equal to the posterior probabilities.
    #' @param scale The scale used for the posterior probability. Either log-scale ("log", default) or natural probabilities summing up to 1 ("prob").
    #' @param prior User-specified prior group probabilities in the new data. If NULL (default), prior probabilities are computed from the learning set.
    #' @param control a list for controlling the optimization. See [PLN()] for details.
    #' @param envir Environment in which the prediction is evaluated
    predict = function(newdata,
                       type = c("posterior", "response", "scores"),
                       scale = c("log", "prob"),
                       prior = NULL,
                       control = PLN_param(backend="nlopt"), envir = parent.frame()) {

      type <- match.arg(type)

      if (type == "scores") scale <- "prob"
      scale <- match.arg(scale)

      ## Extract the model matrices from the new data set with initial formula
      args <- extract_model(call = call("PLNLDA", formula = private$formula, data = newdata), envir = envir)
      ## Remove intercept to prevent interference with binary coding of the grouping factor
      args$X <- args$X[ , colnames(args$X) != "(Intercept)", drop = FALSE]

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
        ve_step <- self$optimize_vestep(X, args$O, args$Y, args$w, control = control)
        cond.log.lik[, k] <- ve_step$Ji
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
        rownames(scores) <- rownames(newdata)
        return(scores)
      }

      if (type == "response") {
        res <- apply(res, 1, which.max)
        res <- factor(groups[res], levels = groups)
        names(res) <- rownames(newdata)
      }
      res
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print methods -------------------------
    #' @description User friendly print method
    show = function() {
      super$show(paste0("Linear Discriminant Analysis for Poisson Lognormal distribution\n"))
      cat("* Additional fields for LDA\n")
      cat("    $percent_var, $corr_map, $scores, $group_means\n")
      cat("* Additional S3 methods for LDA\n")
      cat("    plot.PLNLDAfit(), predict.PLNLDAfit()\n")
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Other functions -----------------------
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    C        = NULL,
    P        = NULL,
    Mu       = NULL,
    grouping = NULL,
    svdLDA   = NULL
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field rank the dimension of the current model
    rank = function() {nlevels(private$grouping) - 1},
    #' @field nb_param number of parameters in the current PLN model
    nb_param = function() {self$p * (self$d + self$rank)},
    #' @field model_par a list with the matrices associated with the estimated parameters of the PLN model: B (covariates), Sigma (latent covariance), C (latent loadings), P (latent position) and Mu (group means)
    model_par = function() {
      par <- super$model_par
      par$C  <- private$C
      par$P  <- private$P
      par$Mu <- private$Mu
      par
    },
    #' @field percent_var the percent of variance explained by each axis
    percent_var = function() {
      eigen.val <- private$svdLDA$d[1:self$rank]^2
      setNames(round(eigen.val/sum(eigen.val)*self$R_squared,4), paste0("LD", 1:self$rank))
    },
    #' @field corr_map a matrix of correlations to plot the correlation circles
    corr_map = function() {
      corr <- cor(private$P, self$scores)
      rownames(corr) <- rownames(private$C)
      colnames(corr) <- paste0("LD", 1:self$rank)
      corr
    },
    #' @field scores a matrix of scores to plot the individual factor maps
    scores     = function() {
      scores <- private$P %*% t(t(private$svdLDA$u[, 1:self$rank]) * private$svdLDA$d[1:self$rank])
      rownames(scores) <- rownames(private$M)
      colnames(scores) <- paste0("LD", 1:self$rank)
      scores
    },
    #' @field group_means a matrix of group mean vectors in the latent space.
    group_means = function() {
      self$model_par$Mu
    }
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNLDAfit_diagonal
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a LDA framework with diagonal covariance
#'
#' @description The function [PLNLDA()] produces an instance of an object with class [`PLNLDAfit`].
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by [PLNfit()], the [plot()] method for
#' LDA visualization and [predict()] method for prediction
#'
#' @param grouping a factor specifying the class of each observation used for discriminant analysis.
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "diagonal"))
#' class(myPLNLDA)
#' print(myPLNLDA)
#' }
PLNLDAfit_diagonal <- R6Class(
  classname = "PLNLDAfit_diagonal",
  inherit = PLNLDAfit,
  public  = list(
    #' @description Initialize a [`PLNfit`] model
    initialize = function(grouping, responses, covariates, offsets, weights, formula, control) {
      super$initialize(grouping, responses, covariates, offsets, weights, formula, control)
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize_diagonal, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_diagonal
    }
  ),
  private = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    torch_elbo = function(data, params, index=torch_tensor(1:self$n)) {
      S2 <- torch_square(params$S[index])
      Z <- data$O[index] + params$M[index] + torch_matmul(data$X[index], params$B)
      res <- .5 * sum(data$w[index]) * sum(torch_log(private$torch_sigma_diag(data, params, index))) +
        sum(data$w[index,NULL] * (torch_exp(Z + .5 * S2) - data$Y[index] * Z -  .5 * torch_log(S2)))
      res
    },

    torch_sigma_diag = function(data, params, index=torch_tensor(1:self$n)) {
      torch_sum(data$w[index,NULL] * (torch_square(params$M[index]) + torch_square(params$S[index])), 1) / sum(data$w[index])
    },

    torch_Sigma = function(data, params, index=torch_tensor(1:self$n)) {
      torch_diag(private$torch_sigma_diag(data, params, index))
    },

    torch_vloglik = function(data, params) {
      S2 <- torch_pow(params$S, 2)
      omega_diag <- torch_pow(private$torch_sigma_diag(data, params), -1)

      Ji <- .5 * self$p - rowSums(.logfactorial(as.matrix(data$Y))) + as.numeric(
        .5 * sum(torch_log(omega_diag)) +
          torch_sum(data$Y * params$Z - params$A + .5 * torch_log(S2), dim = 2) -
          .5 * torch_matmul(torch_pow(params$M, 2) + S2, omega_diag)
      )
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    }
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF TORCH METHODS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ),
  active = list(
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"diagonal"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNLDAfit_diagonal
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNLDAfit_spherical
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a LDA framework with spherical covariance
#'
#' @description The function [PLNLDA()] produces an instance of an object with class [`PLNLDAfit`].
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by [PLNfit()], the [plot()] method for
#' LDA visualization and [predict()] method for prediction
#'
#' @param grouping a factor specifying the class of each observation used for discriminant analysis.
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#'
#' @rdname PLNfit_diagonal
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical"))
#' class(myPLNLDA)
#' print(myPLNLDA)
#' }
PLNLDAfit_spherical <- R6Class(
  classname = "PLNLDAfit_spherical",
  inherit = PLNLDAfit,
  public  = list(
    #' @description Initialize a [`PLNfit`] model
    initialize = function(grouping, responses, covariates, offsets, weights, formula, control) {
      super$initialize(grouping, responses, covariates, offsets, weights, formula, control)
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize_spherical, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_diagonal
    }
  ),
  private = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torch_elbo = function(data, params, index=torch_tensor(1:self$n)) {
      S2 <- torch_square(params$S[index])
      Z <- data$O[index] + params$M[index] + torch_mm(data$X[index], params$B)
      res <- .5 * sum(data$w[index]) * self$p * torch_log(private$torch_sigma2(data, params, index)) -
        sum(data$w[index,NULL] * (data$Y[index] * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
      res
    },

    torch_sigma2 = function(data, params, index=torch_tensor(1:self$n)) {
      sum(data$w[index, NULL] * (torch_square(params$M) + torch_square(params$S))) / (sum(data$w) * self$p)
    },

    torch_Sigma = function(data, params, index=torch_tensor(1:self$n)) {
      torch_eye(self$p) * private$torch_sigma2(data, params, index)
    },

    torch_vloglik = function(data, params) {
      S2 <- torch_pow(params$S, 2)
      sigma2 <- private$torch_sigma2(data, params)
      Ji <- .5 * self$p - rowSums(.logfactorial(as.matrix(data$Y))) + as.numeric(
        torch_sum(data$Y * params$Z - params$A + .5 * torch_log(S2/sigma2) - .5 * (torch_pow(params$M, 2) + S2)/sigma2, dim = 2)
      )
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    }
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ),
  active = list(
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"spherical"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNLDAfit_spherical
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
