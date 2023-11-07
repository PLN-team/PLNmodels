## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNblockfit_diagonal
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
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
## Parameters specific to PLNnetwork-fit methods
#' @param blocks blocks indicators for regrouping the variables/responses (dimension of the residual covariance)
#'
#' @rdname PLNblockfit
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLNblock(Abundance ~ 1, data = trichoptera, blocks = 1:5)
#' class(myPLN)
#' print(myPLN)
#' }
PLNblockfit <- R6Class(
  classname = "PLNblockfit",
  inherit = PLNfit,
  public  = list(
    #' @description Initialize a [`PLNblockfit`] model
    initialize = function(blocks, responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)

      ## Initial memberships/blocks
      private$q <- ncol(blocks)

      ## Overwrite PLNfit Variational parameters (dimension q)
      private$M   <- private$M %*% Z
      private$S   <- matrix(.1, self$n, self$q)
      private$Tau <- t(Z)

      ## Setup of the optimization backend
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize_block, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_block
    }
  ),
  private = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    torch_ELBO = function(data, params, Tau, index=torch_tensor(1:self$n)){
      S2 <- torch_square(params$S[index])
      XB <- torch_mm(data$X[index], params$B[index])
      A <- torch_exp(data$O[index] + XB) * torch_mm(torch_exp(params$M[index] + S2 / 2), Tau)

      res <- sum(data$w[index,NULL] * (data$Y[index] * (data$O[index] + XB + torch_mm(params$M[index], Tau))  - A + .5 * torch_log(S2)))
      res
    },

    torch_Tau = function(data, params) {
      A <- torch_t(torch_exp(params$M + .5 * torch_square(params$S)))
      B <- torch_exp(torch_mm(data$X, params$B) + data$O)
      log_Alpha <- log(torch_unsqueeze(torch_mean(Tau, 2), 2))
      with_no_grad({
        rho <- torch_mm(log_Alpha, torch_ones(c(1,self$p))) +
          torch_mm(torch_t(params$M), data$Y) - torch_mm(A, B)
      })
      torch_clamp(torch:::torch_softmax(rho, 1), 1e-3, 1 - 1e-3)
    },

    torch_vloglik = function(data, params) {
      S2 <- torch_square(params$S)
      Ji <- .5 * self$p - rowSums(.logfactorial(as.matrix(data$Y))) + as.numeric(
        .5 * torch_logdet(params$Omega) +
          torch_sum(data$Y * (params$mu + torch_mm(params$M, params$Tau)) - params$A + .5 * torch_log(S2), dim = 2) -
          .5 * torch_sum(torch_mm(params$M, params$Omega) * params$M + S2 * torch_diag(params$Omega), dim = 2) -
            torch_sum(torch:::torch_xlogy(params$Tau, params$Tau)) +
            torch_sum(torch_mm(torch_t(torch_log(params$Alpha)), params$Tau))
      )
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    },

    #' @import torch
    torch_optimize = function(data, params, config) {

      ## Conversion of data and parameters to torch tensors (pointers)
      data   <- lapply(data, torch_tensor)                         # list with Y, X, O, w
      Tau <- torch_tensor(params$Tau); params$Tau <- NULL
      params <- lapply(params, torch_tensor, requires_grad = TRUE) # list with B, M, S

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
          loss <- - private$torch_elbo(data, params, Tau, index) # compute current ELBO
          loss$backward()                   # backward propagation
          optimizer$step()                  # optimization
        }

        ## Update Tau only every xx iterations
        if(iterate %% it_update_tau == 0){
          Tau   <- private$torch_Tau(data, params)
          Alpha <- torch_unsqueeze(torch_mean(Tau, 2), 2)
          nUpdatesTauDone = nUpdatesTauDone + 1
          last_iupdate = last_iupdate + 1
        }

        ## assess convergence
        objective[iterate + 1] <- - loss$item()
        B_new   <- as.numeric(optimizer$param_groups[[1]]$params$B)
        delta_f <- abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1])
        delta_x <- sum(abs(B_old - B_new))/sum(abs(B_new))

        ## Error message if objective diverges
        if (!is.finite(loss$item())) {
          stop(sprintf("The ELBO diverged during the optimization procedure.\nConsider using:\n* a different optimizer (current optimizer: %s)\n* a smaller learning rate (current rate: %.3f)\nwith `control = PLN_param(config_optim = list(algorithm = ..., lr = ...))`",
                       config$algorithm, config$lr))
        }

        ## display progress
        if (config$trace >  1 && (iterate %% 50 == 0))
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

      params$Tau   <- private$torch_Tau(data, params)
      params$Alpha <- torch_unsqueeze(torch_mean(params$Tau, 2), 2)
      params$Sigma <- private$torch_Sigma(data, params)
      params$Omega <- private$torch_Omega(data, params)
      params$mu    <- data$O + torch_mm(data$X, params$B)
      params$A     <- torch_exp(params$mu) * torch_mm(torch_exp(params$M + torch_square(params$S) / 2), params$Tau)

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
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + .5 * self$q * (self$q + 1) + self$q - 1)},
    #' @field nb_block number blocks of variables (dimension of the residual covariance)
    nb_block   = function() {as.integer(self$q)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"block"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNblockfit_diagonal
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
