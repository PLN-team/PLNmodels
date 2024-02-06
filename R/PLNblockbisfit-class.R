## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNblockbisfit
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit with block-wise residual covariance
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization.
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables (group scale)
#' @param Omega precision matrix of the latent variables (group scale). Inverse of Sigma.
#' @param Tau matrix of posterior probabilities for block belonging
#'
## Parameters specific to PLNnetwork-fit methods
#' @param blocks blocks indicators for regrouping the variables/responses (dimension of the residual covariance)
#'
#' @rdname PLNblockbisfit
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNblockbis(Abundance ~ 1, data = trichoptera, nb_blocks = 1:15)
#' class(myPLN)
#' print(myPLN)
#' }
PLNblockbisfit <- R6Class(
  classname = "PLNblockbisfit",
  inherit = PLNfit,
  public  = list(
    #' @description Initialize a [`PLNblockbisfit`] model
    initialize = function(blocks, responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)

      ## Initial memberships/blocks
      ## Overwrite PLNfit Variational parameters (dimension q)
      private$Delta   <- private$S
      private$Mu   <- private$M
      private$M   <- private$M %*% blocks
      private$S   <- private$S %*% blocks
      private$Tau <- t(blocks)

      ## Setup of the optimization backend
      private$optimizer$main <-
        switch(control$backend,
               "torch" = private$torch_optimize,
               "nlopt" =  nlopt_optimize_block,
               "nlopt-vem" =  optimize_plnblockbis
        )
      private$optimizer$vestep <- nlopt_optimize_vestep_block
    },
    #' @description Update fields of a [`PLNnetworkfit`] object
    #' @param B matrix of regression matrix
    #' @param Sigma variance-covariance matrix of the latent variables (group scale)
    #' @param Omega precision matrix of the latent variables. Inverse of Sigma.(group scale)
    #' @param D variance-covariance matrix of the latent variables (species scale)
    #' @param M     matrix of mean vectors for the variational approximation (group scale)
    #' @param S     matrix of variance vectors for the variational approximation (group scale)
    #' @param Mu     matrix of mean vectors for the variational approximation (species scale)
    #' @param Delta     matrix of variance vectors for the variational approximation (species scale)
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param Tau matrix of posterior probabilities for block belonging
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2    approximate R^2 goodness-of-fit criterion
    #' @param monitoring a list with optimization monitoring quantities
    update = function(Tau=NA, B=NA, Sigma=NA, Omega=NA, D=NA, M=NA, S=NA, Mu=NA, Delta=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
      super$update(B = B, Sigma = Sigma, Omega = Omega, M=M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
      if (!anyNA(Tau)) private$Tau <- Tau
      private$Mu <- Mu
      private$Delta <- Delta
      private$D <- D
    },

    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   params = list(B = private$B, M = private$M, S = private$S,
                                 Mu = private$Mu, Delta = private$Delta, Tau = private$Tau),
                   config = config)
      optim_out <- do.call(private$optimizer$main, args)

      # browser()
      do.call(self$update, optim_out)
    }

  ),
  private = list(
    D = NA,
    Tau = NA, # variational parameters for the block memberships
    Mu = NA, # variational parameter for species mean
    Delta = NA, # variational parameter for species variance - covariance

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    torch_elbo = function(data, params, index=torch_tensor(1:self$n)){
      S2 <- torch_square(params$S[index])
      XB <- torch_mm(data$X[index], params$B)
      A  <- torch_exp(data$O[index] + XB) * torch_mm(torch_exp(params$M[index] + S2 / 2), private$Tau)
      res <- self$n/2 * torch_logdet(private$torch_Sigma(data, params, index)) +
        torch_sum((A - data$Y[index] * (data$O[index] + XB + torch_mm(params$M[index], private$Tau)))) -
        .5 * torch_sum(torch_log(S2))
      res
    },

    torch_update_Tau = function(data, params) {
      A <- torch_t(torch_exp(params$M + .5 * torch_square(params$S)))
      B <- torch_exp(torch_mm(data$X, params$B) + data$O)
      log_Alpha <- torch_log(torch_unsqueeze(torch_mean(private$Tau, 2), 2))
      with_no_grad({
        rho <- torch_mm(log_Alpha, torch_ones(c(1,self$p))) +
          torch_mm(torch_t(params$M), data$Y) - torch_mm(A, B)
      })
      private$Tau <- torch_clamp(torch:::torch_softmax(rho, 1), 1e-3, 1 - 1e-3)
    },

    torch_vloglik = function(data, params) {
      S2 <- torch_square(params$S)
      log_Alpha <- log(torch_unsqueeze(torch_mean(params$Tau, 2), 2))
      Ji <- .5 * self$p - rowSums(.logfactorial(as.matrix(data$Y))) + as.numeric(
        .5 * torch_logdet(params$Omega) +
          torch_sum(data$Y * params$Z - params$A, dim = 2) +
          .5 * torch_sum(torch_log(S2), dim = 2) -
          .5 * torch_sum(torch_mm(params$M, params$Omega) * params$M + S2 * torch_diag(params$Omega), dim = 2) -
          torch_sum(torch:::torch_xlogy(params$Tau, params$Tau)) +
          torch_sum(torch_mm(torch_t(log_Alpha), params$Tau))
      )
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    },

    #' @import torch
    torch_optimize = function(data, params, config) {

      ## Conversion of data and parameters to torch tensors (pointers)
      data   <- lapply(data, torch_tensor)                         # list with Y, X, O, w and Tau0
      params <- lapply(params, torch_tensor, requires_grad = TRUE) # list with B, M, S
      private$Tau <- data$Tau

      ## Initialize optimizer
      optimizer <- switch(config$algorithm,
                          "RPROP"   = optim_rprop(params  , lr = config$lr, etas = config$etas, step_sizes = config$step_sizes),
                          "RMSPROP" = optim_rmsprop(params, lr = config$lr, weight_decay = config$weight_decay, momentum = config$momentum, centered = config$centered),
                          "ADAM"    = optim_adam(params   , lr = config$lr, weight_decay = config$weight_decay),
                          "ADAGRAD" = optim_adagrad(params, lr = config$lr, weight_decay = config$weight_decay)
      )

      ## Optimization loop
      status <- 5
      num_epoch  <- config$num_epoch
      num_batch  <- config$num_batch
      batch_size <- floor(self$n/num_batch)
      objective <- double(length = config$num_epoch + 1)
      for (iterate in 1:num_epoch) {
        B_old <- as.numeric(optimizer$param_groups[[1]]$params$B)

        # rearrange the data each epoch
        permute <- torch::torch_randperm(self$n) + 1L
        for (batch_idx in 1:num_batch) {
          # here index is a vector of the indices in the batch
          index <- permute[(batch_size*(batch_idx - 1) + 1):(batch_idx*batch_size)]

          ## Optimization
          optimizer$zero_grad() # reinitialize gradients
          loss <- private$torch_elbo(data, params, index) # compute current ELBO
          loss$backward()                   # backward propagation
          optimizer$step()                  # optimization
        }

        ## Update Tau only every xx iterations
        if (iterate %% config$it_update_tau == 0) {
          private$torch_update_Tau(data, params)
        }

        ## assess convergence
        objective[iterate + 1] <- loss$item()
        B_new   <- as.numeric(optimizer$param_groups[[1]]$params$B)
        delta_f <- abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1])
        delta_x <- sum(abs(B_old - B_new))/sum(abs(B_new))

        ## Error message if objective diverges
        if (!is.finite(loss$item())) {
          stop(sprintf("The ELBO diverged during the optimization procedure.\nConsider using:\n* a different optimizer (current optimizer: %s)\n* a smaller learning rate (current rate: %.3f)\nwith `control = PLN_param(config_optim = list(algorithm = ..., lr = ...))`",
                       config$algorithm, config$lr))
        }

        ## display progress
        if (config$trace >=  1 && (iterate %% 50 == 0))
          cat('\niteration: ', iterate, 'objective', objective[iterate + 1],
              'delta_f'  , round(delta_f, 6), 'delta_x', round(delta_x, 6))

        ## Check for convergence
        if (delta_f < config$ftol_rel) status <- 3
        if (delta_x < config$xtol_rel) status <- 4
        if (status %in% c(3,4)) {
          objective <- objective[1:iterate + 1]
          break
        }
      }

      params$Tau   <- private$Tau
      params$Sigma <- private$torch_Sigma(data, params)
      params$Omega <- private$torch_Omega(data, params)
      params$Z     <- data$O + torch_mm(data$X, params$B) + torch_mm(params$M, params$Tau)
      params$A     <- torch_exp(data$O + torch_mm(data$X, params$B)) * torch_mm(torch_exp(params$M + torch_square(params$S) / 2), params$Tau)

      out <- lapply(params, as.matrix)
      out$Ji <- private$torch_vloglik(data, params)
      out$monitoring <- list(
        objective  = objective,
        iterations = iterate,
        status     = status,
        backend = "torch"
      )
      out
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF TORCH METHODS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ),
  active = list(
    ####################################
    #' @field blockbis_model_par to add D as a model par
    blockbis_model_par  = function() {list(D = private$D)},
    #' @field blockbis_var_par to add Mu, Delta as var par
    blockbis_var_par  = function() {list(Mu = private$Mu, Delta = private$Delta)},
    ####################################

    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + .5 * self$q * (self$q + 1) + self$q - 1) + self$p},
    #' @field nb_block number blocks of variables (dimension of the residual covariance)
    nb_block   = function() {as.integer(self$q)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"blocks"},
    #' @field loglik (weighted) variational lower bound of the loglikelihood
    loglik     = function() {
      ## this term should not be sumeed over sites/individuals
      pi <- matrix(self$groupProportion, self$q, self$p)
      tau_log_alpha_over_tau <- sum(.xlogy(private$Tau,pi)) - sum(.xlogx(private$Tau))
      sum(self$weights * private$Ji - tau_log_alpha_over_tau) +  tau_log_alpha_over_tau
    },
    #' @field loglik_vec element-wise variational lower bound of the loglikelihood
    loglik_vec = function() {private$Ji},
    #' @field entropy_blocks Entropy of the variational distribution of the block (multinomial)
    entropy_blocks = function() {-sum(.xlogx(private$tau))},
    #' @field entropy_latent Entropy of the variational distribution of the latent vector (Gaussian)
    entropy_latent = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(self$var_par$S2)))},
    #' @field entropy Full entropy of the variational distribution (latent vector + block)
    entropy = function() {self$entropy_latent + self$entropy_blocks},
    #' @field membership block membership for variables
    membership = function() {apply(private$Tau, 2, which.max)},
    #' @field posteriorProb matrix of posterior probabilities for block belonging
    posteriorProb = function() {private$Tau},
    #' @field groupProportion vector of cluster proportions
    groupProportion  = function() {rowMeans(private$Tau)}
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNblockbisfit
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNblockbisfit_sparse
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit with sparse block-wise residual covariance
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization.
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#' @param Tau matrix of posterior probabilities for block belonging
#'
## Parameters specific to PLNblockbis-fit methods
#' @param sparsity tuning parameter for controlling the sparsity level of Omega/Sigma
#' @param blocks blocks indicators for regrouping the variables/responses (dimension of the residual covariance)
#'
#' @rdname PLNblockbisfit_sparse
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNblockbis(Abundance ~ 1, data = trichoptera, nb_blocks = 1:5, sparsity = 0.1)
#' class(myPLN)
#' print(myPLN)
#' }
PLNblockbisfit_sparse <- R6Class(
  classname = "PLNblockbisfit_sparse",
  inherit = PLNblockbisfit,
  private = list(lambda = NA, rho = NA),
  public  = list(
    #' @description Initialize a [`PLNblockbisfit_sparse`] model
    initialize = function(blocks, sparsity, responses, covariates, offsets, weights, formula, control) {
      super$initialize(blocks, responses, covariates, offsets, weights, formula, control)
      private$optimizer$main <- nlopt_optimize_block_sparse
      private$lambda <- sparsity
      private$rho    <- matrix(1, self$q, self$q)
    },
    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   params = list(B = private$B, M = private$M, S = private$S,  Mu = private$Mu, Delta = private$Delta, Tau = private$Tau, rho = private$lambda * private$rho),
                   config = config)
      optim_out <- do.call(private$optimizer$main, args)
      do.call(self$update, optim_out)
    },
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @importFrom Matrix Matrix
    #' @return a square matrix of size `PLNnetworkfit$n`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omega != 0 & !diag(TRUE, ncol(private$Omega))),
        "precision"   = private$Omega,
        "partial_cor" = {
          tmp <- -private$Omega / tcrossprod(sqrt(diag(private$Omega))); diag(tmp) <- 1
          tmp
        }
      )
      ## Enforce sparse Matrix encoding to avoid downstream problems with igraph::graph_from_adjacency_matrix
      ## as it fails when given dsyMatrix objects
      Matrix(net, sparse = TRUE)
    },
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @description plot the latent network.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = layout_in_circle,
                            plot = TRUE) {

      type <- match.arg(type)
      output <- match.arg(output)

      net <- self$latent_network(type)

      if (output == "igraph") {

        G <-  graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- colnames(net)
        }
        ## Nice nodes
        V.deg <- degree(G)/sum(degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- V.deg * 100
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- NA
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- delete.vertices(G, which(degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
    }
  ),
  active = list(
    #' @field penalty tuning parameter for controlling the sparsity level of Omega/Sigma
    penalty         = function() {private$lambda},
    #' @field penalty_weights a matrix of weights controlling the amount of penalty element-wise.
    penalty_weights = function() {private$rho},
    #' @field n_edges number of edges if the network (non null coefficient of the sparse precision matrix)
    n_edges         = function() {sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0)},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + self$n_edges + self$q - 1)},
    #' @field density proportion of non-null edges in the network
    density   = function() {mean(self$latent_network("support"))},
    #' @field pen_loglik variational lower bound of the l1-penalized loglikelihood
    pen_loglik      = function() {self$loglik - private$lambda * sum(abs(private$Omega))},
    #' @field EBIC variational lower bound of the EBIC
    EBIC      = function() {self$BIC - .5 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$p*(self$p - 1)/self$n_edges), 0)},
    #' @field criteria a vector with loglik, penalized loglik, BIC, EBIC, ICL, number of parameters, number of edges and graph density
    criteria  = function() {data.frame(super$criteria, n_edges = self$n_edges, EBIC = self$EBIC, pen_loglik = self$pen_loglik, density = self$density)}
  )
)
