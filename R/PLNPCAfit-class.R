#' An R6 Class to represent a PLNfit in a PCA framework
#'
#' @description The function [PLNPCA()] produces a collection of models which are instances of object with class [`PLNPCAfit`].
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for the methods inherited by  [`PLNfit`] and the [plot()] methods for PCA visualization
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in [`PLNfamily`]
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
## Parameters specific to PLNPCA-fit methods
#' @param rank rank of the PCA (or equivalently, dimension of the latent space)
#'
## Parameters common to many graphical methods
#' @param map the type of output for the PCA visualization: either "individual", "variable" or "both". Default is "both".
#' @param nb_axes scalar: the number of axes to be considered when map = "both". The default is min(3,rank).
#' @param axes numeric, the axes to use for the plot when map = "individual" or "variable". Default it c(1,min(rank))
#' @param ind_cols a character, factor or numeric to define the color associated with the individuals. By default, all variables receive the default color of the current palette.
#' @param var_cols a character, factor or numeric to define the color associated with the variables. By default, all variables receive the default color of the current palette.
#' @param plot logical. Should the plot be displayed or sent back as ggplot object
#' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
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
    ## PRIVATE MEMBERS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      C     = NULL,
      svdCM = NULL,

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## PRIVATE TORCH METHODS FOR RANK-CONSTRAINED OPTIMIZATION
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      torch_elbo_rank_core = function(data, M, S, B, C, index) {
        S2 <- torch_square(S[index])                        # (batch, q)
        C2 <- torch_square(C)                              # (p, q)
        Z  <- data$O[index] +
              torch_mm(M[index], torch_t(C)) +
              torch_mm(data$X[index], B)                   # (batch, p)
        A  <- torch_exp(Z + 0.5 * torch_mm(S2, torch_t(C2)))
        lik_part <- torch_sum(data$w[index, NULL] * (A - data$Y[index] * Z))
        kl_part  <- 0.5 * torch_sum(data$w[index, NULL] *
                      (torch_square(M[index]) + S2 - torch_log(S2) - 1))
        lik_part + kl_part
      },

      torch_elbo_rank = function(data, params, index = torch_tensor(1:self$n)) {
        private$torch_elbo_rank_core(data, params$M, params$S, params$B, params$C, index)
      },

      torch_vloglik_rank = function(data, params) {
        S2 <- torch_square(params$S)
        C2 <- torch_square(params$C)
        Z  <- data$O + torch_mm(params$M, torch_t(params$C)) + torch_mm(data$X, params$B)
        A  <- torch_exp(Z + 0.5 * torch_mm(S2, torch_t(C2)))
        Ji <- - torch_sum(.logfactorial_torch(data$Y), dim = 2) +
              torch_sum(data$Y * Z - A, dim = 2) -
              0.5 * torch_sum(torch_square(params$M) + S2 - torch_log(S2) - 1, dim = 2)
        Ji <- .5 * self$p + as.numeric(Ji$cpu())
        attr(Ji, "weights") <- as.numeric(data$w$cpu())
        Ji
      },

      torch_optimize_rank_core = function(data, params, config, n_obs, loss_fn) {
        optimizer <- switch(config$algorithm,
          "RPROP"   = optim_rprop(params,   lr = config$lr, etas = config$etas, step_sizes = config$step_sizes),
          "RMSPROP" = optim_rmsprop(params, lr = config$lr, weight_decay = config$weight_decay, momentum = config$momentum, centered = config$centered),
          "ADAM"    = optim_adam(params,    lr = config$lr, weight_decay = config$weight_decay),
          "ADAGRAD" = optim_adagrad(params, lr = config$lr, weight_decay = config$weight_decay)
        )

        status     <- 5
        num_epoch  <- config$num_epoch
        num_batch  <- config$num_batch
        batch_size <- floor(n_obs / num_batch)

        objective <- double(length = config$num_epoch + 1)
        for (iterate in 1:num_epoch) {
          permute <- torch::torch_tensor(sample.int(n_obs), dtype = torch_long(), device = config$device)
          for (batch_idx in 1:num_batch) {
            index <- permute[(batch_size * (batch_idx - 1) + 1):(batch_idx * batch_size)]
            optimizer$zero_grad()
            loss <- loss_fn(index)
            loss$backward()
            optimizer$step()
          }

          objective[iterate + 1] <- loss$item()
          delta_f <- abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1])

          if (!is.finite(loss$item())) {
            stop(sprintf(
              "The ELBO diverged during the optimization procedure.\nConsider using:\n* a different optimizer (current optimizer: %s)\n* a smaller learning rate (current rate: %.3f)\nwith `control = PLNPCA_param(backend = 'torch', config_optim = list(algorithm = ..., lr = ...))`",
              config$algorithm, config$lr
            ))
          }

          if (config$trace > 1 && (iterate %% 50 == 1))
            cat('\niteration:', iterate, 'objective', objective[iterate + 1],
                'delta_f', round(delta_f, 6))

          if (delta_f < config$ftol_rel) status <- 3
          if (status %in% c(3, 4)) {
            objective <- objective[seq_len(iterate + 1)]
            break
          }
        }

        list(
          params     = params,
          objective  = objective,
          iterations = iterate,
          status     = status
        )
      },

      torch_optimize_vestep_rank = function(data, params, B, C, config) {
        if (config$trace > 1)
          message(paste("optimizing with device:", config$device))

        n <- nrow(data$Y)
        data   <- lapply(data,   torch_tensor, dtype = torch_float32(), device = config$device)
        params <- lapply(params, torch_tensor, dtype = torch_float32(), requires_grad = TRUE, device = config$device)
        B      <- torch_tensor(B, dtype = torch_float32(), device = config$device)
        C      <- torch_tensor(C, dtype = torch_float32(), device = config$device)

        optim_out <- private$torch_optimize_rank_core(
          data   = data,
          params = params,
          config = config,
          n_obs  = n,
          loss_fn = function(index) {
            private$torch_elbo_rank_core(data, params$M, params$S, B, C, index)
          }
        )
        params_r <- lapply(optim_out$params, function(x) as.matrix(x$cpu()))
        Ji_r <- private$torch_vloglik_rank(data, c(optim_out$params, list(B = B, C = C)))

        list(
          M          = params_r$M,
          S          = params_r$S,
          Ji         = Ji_r,
          monitoring = list(
            objective  = optim_out$objective,
            iterations = optim_out$iterations,
            status     = optim_out$status,
            backend    = "torch"
          )
        )
      },

      torch_optimize_rank = function(data, params, config) {
        if (config$trace > 1)
          message(paste("optimizing with device:", config$device))

        data   <- lapply(data,   torch_tensor, dtype = torch_float32(), device = config$device)
        params <- lapply(params, torch_tensor, dtype = torch_float32(), requires_grad = TRUE, device = config$device)

        optim_out <- private$torch_optimize_rank_core(
          data   = data,
          params = params,
          config = config,
          n_obs  = self$n,
          loss_fn = function(index) {
            private$torch_elbo_rank(data, params, index)
          }
        )

        ## Compute derived quantities on CPU
        params_r <- lapply(optim_out$params, function(x) as.matrix(x$cpu()))
        data_r   <- lapply(data,   function(x) as.matrix(x$cpu()))

        q     <- ncol(params_r$M)
        S2_r  <- params_r$S^2
        C2_r  <- params_r$C^2
        Z_r   <- data_r$O + params_r$M %*% t(params_r$C) + data_r$X %*% params_r$B
        A_r   <- exp(Z_r + 0.5 * S2_r %*% t(C2_r))
        w_r   <- as.numeric(data_r$w)

        wM      <- params_r$M * sqrt(w_r)
        inner_q <- (crossprod(wM) + diag(colSums(S2_r * w_r), nrow = q)) / sum(w_r)
        Sigma_r <- params_r$C %*% inner_q %*% t(params_r$C)
        Omega_r <- params_r$C %*% solve(inner_q) %*% t(params_r$C)

        Ji_r <- .5 * self$p - rowSums(.logfactorial(as.matrix(data_r$Y))) +
                rowSums(data_r$Y * Z_r - A_r) -
                0.5 * rowSums(params_r$M^2 + S2_r - log(S2_r) - 1)
        attr(Ji_r, "weights") <- w_r

        list(
          B          = params_r$B,
          C          = params_r$C,
          M          = params_r$M,
          S          = params_r$S,
          Z          = Z_r,
          A          = A_r,
          Sigma      = Sigma_r,
          Omega      = Omega_r,
          Ji         = Ji_r,
          monitoring = list(
            objective  = optim_out$objective,
            iterations = optim_out$iterations,
            status     = optim_out$status,
            backend    = "torch"
          )
        )
      }
    ),
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public  = list(
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Creation functions ----------------
      #' @description Initialize a [`PLNPCAfit`] object
      initialize = function(rank, responses, covariates, offsets, weights, formula, control) {
        super$initialize(responses, covariates, offsets, weights, formula, control)
        if (control$backend == "torch") {
          private$optimizer$main <- private$torch_optimize_rank
        } else {
          private$optimizer$main <- nlopt_optimize_rank
        }
        private$optimizer$vestep <- nlopt_optimize_vestep_rank
        if (!is.null(control$svdM)) {
          svdM <- control$svdM
        } else {
          svdM <- svd(private$M, nu = rank, nv = self$p)
        }
        ### TODO: check that it is really better than initializing with zeros...
        private$M  <- svdM$u[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank) %*% t(svdM$v[1:rank, 1:rank, drop = FALSE])
        private$S  <- matrix(0.1, self$n, rank)
        private$C  <- svdM$v[, 1:rank, drop = FALSE] %*% diag(svdM$d[1:rank], nrow = rank, ncol = rank)/sqrt(self$n)
      },
      #' @description Update a [`PLNPCAfit`] object
      #' @param M     matrix of mean vectors for the variational approximation
      #' @param C     matrix of PCA loadings (in the latent space)
      #' @param S     matrix of variance vectors for the variational approximation
      #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
      #' @param R2    approximate R^2 goodness-of-fit criterion
      #' @param Z     matrix of latent vectors (includes covariates and offset effects)
      #' @param A     matrix of fitted values
      #' @param monitoring a list with optimization monitoring quantities
      #' @return Update the current [`PLNPCAfit`] object
      update = function(B=NA, Sigma=NA, Omega=NA, C=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
        super$update(B = B, Sigma = Sigma, Omega = Omega, M = M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
        if (!anyNA(C)) private$C <- C
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Optimization ----------------------
      #' @description Call to the C++ optimizer and update of the relevant fields
      optimize = function(responses, covariates, offsets, weights, config) {
        args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                     params = list(B = private$B, C = private$C, M = private$M, S = private$S),
                     config = config)
        optim_out <- do.call(private$optimizer$main, args)
        do.call(self$update, optim_out)
      },

      #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (C, B). Intended to position new data in the latent space for further use with PCA.
      #' @return A list with three components:
      #'  * the matrix `M` of variational means,
      #'  * the matrix `S2` of variational variances
      #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
      optimize_vestep = function(covariates, offsets, responses, weights = rep(1, self$n),
                                 control = PLNPCA_param(backend = "nlopt")) {

        # problem dimension
        n <- nrow(responses); p <- ncol(responses); q <- self$rank

        # turn offset vector to offset matrix
        offsets <- as.matrix(offsets)
        if (ncol(offsets) == 1) offsets <- matrix(offsets, nrow = n, ncol = p)

        ## Not completely naive starting values for M: SVD on the residuals of
        ## a linear regression on the log-counts (+1 to deal with 0s)
        log_responses <- log(responses+1)
        residuals <- lm.wfit(x = covariates, y = log_responses, w = weights, offset = offsets)$residuals
        svd_residuals <- svd(residuals, nu = q, nv = p)
        M_init <- svd_residuals$u[, 1:q, drop = FALSE] %*% diag(svd_residuals$d[1:q], nrow = q, ncol = q) %*% t(svd_residuals$v[1:q, 1:q, drop = FALSE])

        ## Initialize the variational parameters with the appropriate new dimension of the data
        args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                     ## Initialize the variational parameters with the new dimension of the data
                     params = list(M = M_init, S = matrix(.1, n, q)),
                     B = private$B,
                     C = private$C,
                     config = control$config_optim)
        vestep_optimizer <- if (control$backend == "torch") private$torch_optimize_vestep_rank else private$optimizer$vestep
        optim_out <- do.call(vestep_optimizer, args)
        optim_out
      },

      #' @description Project new samples into the PCA space using one VE step
      #' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
      #' @param control a list for controlling the optimization. See [PLN()] for details.
      #' @param envir Environment in which the projection is evaluated
      #' @return
      #'  * the named matrix of scores for the newdata, expressed in the same coordinate system as `self$scores`
      project = function(newdata, control = PLNPCA_param(), envir = parent.frame()) {

        ## Extract the model matrices from the new data set with initial formula
        args <- extract_model(call("PLNPCA", formula = private$formula, data = newdata), envir)

        ## Compute latent positions of the new samples
        M <- self$optimize_vestep(covariates = args$X, offsets = args$O, responses = args$Y,
                            weights = args$w, control = control)$M
        latent_pos <- t(tcrossprod(self$model_par$C, M)) %>% scale(center = TRUE, scale = FALSE)

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
      setVisualization = function(scale.unit = FALSE) {
        private$svdCM <- svd(scale(self$latent_pos, TRUE, scale.unit), nv = self$rank)
      },

      #' @description Update R2, fisher, std_err fields and set up visualization
      #' @param config_optim a list for controlling the optimizer (either "nlopt" or "torch" backend). See details
      #' @param config_post a list for controlling the post-treatments (optional bootstrap, jackknife, R2, etc.). See details
      #' @details The list of parameters `config_post` controls the post-treatment processing, with the following entries:
      #' * jackknife boolean indicating whether jackknife should be performed to evaluate bias and variance of the model parameters. Default is FALSE.
      #' * bootstrap integer indicating the number of bootstrap resamples generated to evaluate the variance of the model parameters. Default is 0 (inactivated).
      #' * variational_var boolean indicating whether variational Fisher information matrix should be computed to estimate the variance of the model parameters (highly underestimated). Default is FALSE.
      #' * rsquared boolean indicating whether approximation of R2 based on deviance should be computed. Default is TRUE
      #' * trace integer for verbosity. should be > 1 to see output in post-treatments
      postTreatment = function(responses, covariates, offsets, weights, config_post, config_optim, nullModel) {
        super$postTreatment(responses, covariates, offsets, weights, config_post, config_optim, nullModel)
        colnames(private$C) <- colnames(private$M) <- 1:self$q
        rownames(private$C) <- colnames(responses)
        self$setVisualization()
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------

      #' @description Plot the factorial map of the PCA
      # @inheritParams plot.PLNPCAfit
      #' @param cols a character, factor or numeric to define the color associated with the individuals. By default, all individuals receive the default color of the current palette.
      #' @return a [`ggplot2::ggplot`] graphic
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
      #' @return a [`ggplot2::ggplot`] graphic
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
    ##  ACTIVE BINDINGS ----
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    active = list(
      #' @field rank the dimension of the current model
      rank = function() {self$q},
      #' @field vcov_model character: the model used for the residual covariance
      vcov_model = function() {"rank"},
      #' @field nb_param number of parameters in the current PLN model
      nb_param = function() {self$p * (self$d + self$q) - self$q * (self$q - 1)/2},
      #' @field entropy entropy of the variational distribution
      entropy  = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(self$var_par$S2)))},
      #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
      latent_pos = function() {tcrossprod(private$M, private$C)},
      #' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: B (covariates), Sigma (covariance), Omega (precision) and C (loadings)
      model_par = function() {
        par <- super$model_par
        par$C <- private$C
        par
      },
      #' @field percent_var the percent of variance explained by each axis
      percent_var = function() {
        eigen.val <- private$svdCM$d[1:self$rank]^2
        setNames(round(eigen.val/sum(eigen.val)*self$R_squared,4), paste0("PC", 1:self$rank))
      },
      #' @field corr_circle a matrix of correlations to plot the correlation circles
      corr_circle = function() {
        corr <- private$svdCM$v[, 1:self$rank] * matrix(private$svdCM$d[1:self$rank], byrow = TRUE, nrow = self$p, ncol = self$rank)
        corr <- corr/sqrt(rowSums(corr^2))
        rownames(corr) <- rownames(private$Sigma)
        colnames(corr) <- paste0("PC", 1:self$rank)
        corr
      },
      #' @field scores a matrix of scores to plot the individual factor maps (a.k.a. principal components)
      scores     = function() {
        scores <- private$svdCM$u[, 1:self$rank] * matrix(private$svdCM$d[1:self$rank], byrow = TRUE, nrow = self$n, ncol = self$rank)
        rownames(scores) <- rownames(private$M)
        colnames(scores) <- paste0("PC", 1:self$rank)
        scores
      },
      #' @field rotation a matrix of rotation of the latent space
      rotation   = function() {
        rotation <- private$svdCM$v[, 1:self$rank, drop = FALSE]
        rownames(rotation) <- rownames(private$Sigma)
        colnames(rotation) <- paste0("PC", 1:self$rank)
        rotation
      },
      #' @field eig description of the eigenvalues, similar to percent_var but for use with external methods
      eig = function() {
        eigen.val <- private$svdCM$d[1:self$rank]^2
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
        coord  <- private$svdCM$v[, 1:self$rank] * matrix(private$svdCM$d[1:self$rank], ncol = self$rank, nrow = self$p, byrow = TRUE)
        ## coord[j, k] = d[k] * v[j, k]
        var_sd <- sqrt(rowSums(coord^2))
        coord  <- coord / var_sd
        cor    <- coord
        cos2 <- cor^2
        contrib <- 100 * private$svdCM$v[, 1:self$rank, drop = FALSE]^2
        dimnames(coord) <- dimnames(cor) <- dimnames(cos2) <- dimnames(contrib) <- list(rownames(private$Sigma), paste0("Dim.", 1:self$rank))
        list(coord   = coord,
             cor     = cor,
             cos2    = cos2,
             contrib = contrib)
      },
      #' @field ind a list of data frames with PCA results for the individuals: `coord` (coordinates of the individuals), `cos2` (Cosine of the individuals), `contrib` (contributions of individuals to an axis inertia) and `dist` (distance of individuals to the origin).
      ind = function() {
        coord  <- private$svdCM$u[, 1:self$rank] * matrix(private$svdCM$d[1:self$rank], ncol = self$rank, nrow = self$n, byrow = TRUE)
        ## coord[i, k] = d[k] * v[i, k]
        dist_origin <- sqrt(rowSums(coord^2))
        cos2 <- coord^2 / dist_origin^2
        contrib <- 100 * private$svdCM$u[, 1:self$rank, drop = FALSE]^2
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
