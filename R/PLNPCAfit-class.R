#' An R6 Class to represent a PLNfit in a PCA framework
#'
#' @description The function [PLNPCA()] produces a collection of models which are instances of object with class [`PLNPCAfit`].
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by  [`PLNfit`] and the [plot()] methods for PCA visualization
#'
## Parameters common to many PLNPCAfit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param rank rank of the PCA (or equivalently, dimension of the latent space)
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call and used for predictions.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param type approximation scheme to compute the fisher information matrix. Either `wald` (default) or `louis`. `type = "louis"` results in smaller confidence intervals.
## Parameters common to many PLNLDAfit graphical methods
#' @param map the type of output for the PCA visualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param ind_cols a character, factor or numeric to define the color associated with the individuals. By default, all variables receive the default color of the current palette.
#' @param var_cols a character, factor or numeric to define the color associated with the variables. By default, all variables receive the default color of the current palette.
#' @param plot logical. Should the plot be displayed or sent back as ggplot object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
#'
#'
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
PLNPCAfit <- R6Class(
    classname = "PLNPCAfit",
    inherit = PLNfit,
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public  = list(
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Creation functions ----------------
      #' @description Initialize a [`PLNPCAfit`] object
      initialize = function(rank, responses, covariates, offsets, weights, formula, xlevels, control) {
        super$initialize(responses, covariates, offsets, weights, formula, xlevels, control)
        if (!is.null(control$svdM)) {
          svdM <- control$svdM
        } else {
          svdM <- svd(private$M, nu = rank, nv = self$p)
        }
        private$M  <- svdM$u[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank) %*% t(svdM$v[1:rank, 1:rank, drop = FALSE])
        private$S2 <- matrix(0.1, self$n, rank)
        private$B  <- svdM$v[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank)/sqrt(self$n)
        private$covariance <- "rank"
      },
      #' @description Update a [`PLNPCAfit`] object
      #' @param Theta matrix of regression matrix
      #' @param Sigma variance-covariance matrix of the latent variables
      #' @param M     matrix of mean vectors for the variational approximation
      #' @param B     matrix of PCA loadings (in the latent space)
      #' @param S2    matrix of variance vectors for the variational approximation
      #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
      #' @param R2    approximate R^2 goodness-of-fit criterion
      #' @param Z     matrix of latent vectors (includes covariates and offset effects)
      #' @param A     matrix of fitted values
      #' @param monitoring a list with optimization monitoring quantities
      #' @return Update the current [`PLNPCAfit`] object
      update = function(Theta=NA, Sigma=NA, B=NA, M=NA, S2=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
        super$update(Theta = Theta, Sigma = Sigma, M = M, S2 = S2, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
        if (!anyNA(B)) private$B <- B
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Optimization ----------------------
      #' @description Call to the C++ optimizer and update of the relevant fields
      optimize = function(responses, covariates, offsets, weights, control) {
        ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
        opts <- control
        opts$xtol_abs <- list(Theta = 0, B = 0, M = 0, S = control$xtol_abs)
        optim_out <- cpp_optimize_rank(
          list(
            Theta = private$Theta,
            B = private$B,
            M = private$M,
            S = sqrt(private$S2)
          ),
          responses, covariates, offsets, weights, opts
        )

        Ji <- optim_out$loglik
        attr(Ji, "weights") <- weights
        self$update(
          B     = optim_out$B,
          Theta = optim_out$Theta,
          Sigma = optim_out$Sigma,
          M     = optim_out$M,
          S2    = (optim_out$S)**2,
          A     = optim_out$A,
          Z     = optim_out$Z,
          Ji    = Ji,
          monitoring = list(
            iterations = optim_out$iterations,
            status     = optim_out$status,
            message    = statusToMessage(optim_out$status))
        )
      },

      #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (B, Theta). Intended to position new data in the latent space for further use with PCA.
      #' @return A list with three components:
      #'  * the matrix `M` of variational means,
      #'  * the matrix `S2` of variational variances
      #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
      VEstep = function(covariates, offsets, responses, weights = rep(1, nrow(responses)), control = list()) {

        # problem dimension
        n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates); q <- self$rank

        # turn offset vector to offset matrix
        offsets <- as.matrix(offsets)
        if (ncol(offsets) == 1) offsets <- matrix(offsets, nrow = n, ncol = p)

        ## define default control parameters for optim and overwrite by user defined parameters
        control$covariance <- self$vcov_model
        control <- PLN_param(control, n, p)

        VEstep_optimizer  <- cpp_optimize_vestep_rank

        ## Not completely naive starting values for M: SVD on the residuals of
        ## a linear regression on the log-counts (+1 to deal with 0s)
        log_responses <- log(responses+1)
        residuals <- lm.wfit(x = covariates, y = log_responses, w = weights, offset = offsets)$residuals
        svd_residuals <- svd(residuals, nu = q, nv = p)
        M_init <- svd_residuals$u[, 1:q, drop = FALSE] %*% diag(svd_residuals$d[1:q], nrow = q, ncol = q) %*% t(svd_residuals$v[1:q, 1:q, drop = FALSE])

        ## Initialize the variational parameters with the appropriate new dimension of the data
        optim_out <- VEstep_optimizer(
          list(
            ## M = matrix(0, n, q),
            M = M_init,
            S = matrix(sqrt(0.1), n, q)
          ),
          responses,
          covariates,
          offsets,
          weights,
          Theta = self$model_par$Theta,
          B = self$model_par$B,
          control
        )

        Ji <- optim_out$loglik
        attr(Ji, "weights") <- weights

        ## output
        list(M       = optim_out$M,
             S2      = (optim_out$S)**2,
             log.lik = setNames(Ji, rownames(responses)))
      },

      #' @description Project new samples into the PCA space using one VE step
      #' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
      #' @param control a list for controlling the optimization. See [PLN()] for details.
      #' @param envir Environment in which the projection is evaluated
      #' @return
      #'  * the named matrix of scores for the newdata, expressed in the same coordinate system as `self$scores`
      project = function(newdata, control = list(), envir = parent.frame()) {

        ## Extract the model matrices from the new data set with initial formula
        args <- extract_model(call("PLNPCA", formula = private$formula, data = newdata, xlev = private$xlevels), envir)

        ## Compute latent positions of the new samples
        M <- self$VEstep(covariates = args$X, offsets = args$O, responses = args$Y,
                            weights = args$w, control = control)$M
        latent_pos <- t(tcrossprod(self$model_par$B, M)) %>% scale(center = TRUE, scale = FALSE)

        ## Compute scores in the PCA coordinate systems
        scores <- latent_pos %*% self$rotation
        dimnames(scores) <- list(rownames(newdata), paste0("PC", 1:ncol(scores)))

        ## Output
        scores
      },


      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Post treatment --------------------
      #' @description Compute PCA scores in the latent space and update corresponding fields.
      #' @param scale.unit Logical. Should PCA scores be rescaled to have unit variance
      setVisualization = function(scale.unit=FALSE) {
        private$svdBM <- svd(scale(self$latent_pos,TRUE, scale.unit), nv = self$rank)
      },

      #' @description Update R2, fisher, std_err fields and set up visualization
      #' after optimization
      postTreatment = function(responses, covariates, offsets, weights, nullModel) {
        super$postTreatment(responses, covariates, offsets, weights, type = "none", nullModel = nullModel)
        colnames(private$B) <- colnames(private$M) <- 1:self$q
        rownames(private$B) <- colnames(responses)
        if (private$covariance != "spherical") colnames(private$S2) <- 1:self$q
        self$setVisualization()
      },

      #' @description Safely compute the fisher information matrix (FIM)
      #' @param X design matrix used to compute the FIM
      #' @return a sparse matrix with sensible dimension names
      compute_fisher = function(type = c("wald", "louis"), X = NULL) {
        type = match.arg(type)
        if (type == "louis") {
          stop("Louis approximation scheme not available yet for object of class PLNPLCA, use type = \"wald\" instead.")
        }
        super$compute_fisher(type = "wald", X = X)
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------

      #' @description Plot the factorial map of the PCA
      # @inheritParams plot.PLNPCAfit
      #' @param cols a character, factor or numeric to define the color associated with the individuals. By default, all individuals receive the default color of the current palette.
      #' @return a [`ggplot`] graphic
      plot_individual_map = function(axes=1:min(2,self$rank), main = "Individual Factor Map", plot = TRUE, cols = "default") {

        .scores <- data.frame(self$scores[,axes, drop = FALSE])
        colnames(.scores) <- paste("a",1:ncol(.scores),sep = "")
        .scores$labels <- cols
        .scores$names <- rownames(private$M)
        axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

        p <- get_ggplot_ind_map(.scores, axes_label, main)
        if (plot) print(p)
        invisible(p)
      },
      #' @description Plot the correlation circle of a specified axis for a [`PLNLDAfit`] object
      # @inheritParams plot.PLNPCAfit
      #' @param cols a character, factor or numeric to define the color associated with the variables. By default, all variables receive the default color of the current palette.
      #' @return a [`ggplot`] graphic
      plot_correlation_circle = function(axes=1:min(2,self$rank), main="Variable Factor Map", cols = "default", plot=TRUE) {

        ## data frame with correlations between variables and PCs
        correlations <- as.data.frame(self$corr_circle[, axes, drop = FALSE])
        colnames(correlations) <- paste0("axe", 1:length(axes))
        correlations$labels <- cols
        correlations$names  <- rownames(correlations)
        axes_label <- paste(paste("axis",axes), paste0("(", round(100*self$percent_var,3)[axes], "%)"))

        p <- get_ggplot_corr_circle(correlations, axes_label, main)

        if (plot) print(p)
        invisible(p)
      },

      #' @description Plot a summary of the [`PLNPCAfit`] object
      # @inheritParams plot.PLNPCAfit
      #' @importFrom gridExtra grid.arrange arrangeGrob
      #' @importFrom grid nullGrob textGrob
      #' @return a [`grob`] object
      plot_PCA = function(nb_axes = min(3, self$rank), ind_cols = "ind_cols", var_cols = "var_cols", plot = TRUE) {

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
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Print methods ---------------------

      #' @description User friendly print method
      show = function() {
        super$show(paste0("Poisson Lognormal with rank constrained for PCA (rank = ",self$rank,")\n"))
        cat("* Additional fields for PCA\n")
        cat("    $percent_var, $corr_circle, $scores, $rotation, $eig, $var, $ind\n")
        cat("* Additional S3 methods for PCA\n")
        cat("    plot.PLNPCAfit()\n")
      }

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Other methods ---------------------
    ),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE MEMBERS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      B     = NULL,
      svdBM = NULL
    ),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##  ACTIVE BINDINGS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    active = list(
      #' @field rank the dimension of the current model
      rank = function() {self$q},
      #' @field nb_param number of parameters in the current PLN model
      nb_param = function() {self$p * (self$d + self$q) - self$q * (self$q - 1)/2},
      #' @field entropy entropy of the variational distribution
      entropy  = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(private$S2)))},
      #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
      latent_pos = function() {t(tcrossprod(private$B, private$M))},
      #' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: Theta (covariates), Sigma (latent covariance) and B (latent loadings)
      model_par = function() {
        par <- super$model_par
        par$B <- private$B
        par
      },
      #' @field percent_var the percent of variance explained by each axis
      percent_var = function() {
        eigen.val <- private$svdBM$d[1:self$rank]^2
        setNames(round(eigen.val/sum(eigen.val)*self$R_squared,4), paste0("PC", 1:self$rank))
      },
      #' @field corr_circle a matrix of correlations to plot the correlation circles
      corr_circle = function() {
        corr <- private$svdBM$v[, 1:self$rank] * matrix(private$svdBM$d[1:self$rank], byrow = TRUE, nrow = self$p, ncol = self$rank)
        corr <- corr/sqrt(rowSums(corr^2))
        rownames(corr) <- rownames(private$Sigma)
        colnames(corr) <- paste0("PC", 1:self$rank)
        corr
      },
      #' @field scores a matrix of scores to plot the individual factor maps (a.k.a. principal components)
      scores     = function() {
        scores <- private$svdBM$u[, 1:self$rank] * matrix(private$svdBM$d[1:self$rank], byrow = TRUE, nrow = self$n, ncol = self$rank)
        rownames(scores) <- rownames(private$M)
        colnames(scores) <- paste0("PC", 1:self$rank)
        scores
      },
      #' @field rotation a matrix of rotation of the latent space
      rotation   = function() {
        rotation <- private$svdBM$v[, 1:self$rank, drop = FALSE]
        rownames(rotation) <- rownames(private$Sigma)
        colnames(rotation) <- paste0("PC", 1:self$rank)
        rotation
      },
      #' @field eig description of the eigenvalues, similar to percent_var but for use with external methods
      eig = function() {
        eigen.val <- private$svdBM$d[1:self$rank]^2
        matrix(
          c(eigen.val,                                                # eigenvalues
            100 * self$R_squared * eigen.val / sum(eigen.val),        # percentage of variance
            100 * self$R_squared * cumsum(eigen.val) / sum(eigen.val) # cumulative percentage of variance
          ),
          ncol = 3,
          dimnames = list(paste("comp", 1:self$rank), c("eigenvalue", "percentage of variance", "cumulative percentage of variance"))
        )
      },
      #' @field var a list of data frames with PCA results for the variables: `coord` (coordinates of the variables), `cor` (correlation between variables and dimensions), `cos2` (Cosine of the variables) and `contrib` (contributions of the variable to the axes)
      var = function() {
        coord  <- private$svdBM$v[, 1:self$rank] * matrix(private$svdBM$d[1:self$rank], ncol = self$rank, nrow = self$p, byrow = TRUE)
        ## coord[j, k] = d[k] * v[j, k]
        var_sd <- sqrt(rowSums(coord^2))
        coord  <- coord / var_sd
        cor    <- coord
        cos2 <- cor^2
        contrib <- 100 * private$svdBM$v[, 1:self$rank, drop = FALSE]^2
        dimnames(coord) <- dimnames(cor) <- dimnames(cos2) <- dimnames(contrib) <- list(rownames(private$Sigma), paste0("Dim.", 1:self$rank))
        list(coord   = coord,
             cor     = cor,
             cos2    = cos2,
             contrib = contrib)
      },
      #' @field ind a list of data frames with PCA results for the individuals: `coord` (coordinates of the individuals), `cos2` (Cosine of the individuals), `contrib` (contributions of individuals to an axis inertia) and `dist` (distance of individuals to the origin).
      ind = function() {
        coord  <- private$svdBM$u[, 1:self$rank] * matrix(private$svdBM$d[1:self$rank], ncol = self$rank, nrow = self$n, byrow = TRUE)
        ## coord[i, k] = d[k] * v[i, k]
        dist_origin <- sqrt(rowSums(coord^2))
        cos2 <- coord^2 / dist_origin^2
        contrib <- 100 * private$svdBM$u[, 1:self$rank, drop = FALSE]^2
        dimnames(coord) <- dimnames(cos2) <- dimnames(contrib) <- list(rownames(private$M), paste0("Dim.", 1:self$rank))
        names(dist_origin) <- rownames(private$M)
        list(coord   = coord,
             cos2    = cos2,
             contrib = contrib,
             dist    = dist_origin)
      },
      #' @field call Hacky binding for compatibility with factoextra functions
      call = function() {
        list(scale.unit = FALSE)
      }
    )
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##  END OF CLASS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
