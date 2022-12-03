#' An R6 Class to represent a PLNfit in a standard, general framework
#'
#' @description The function [PLN()] fit a model which is an instance of a object with class [`PLNfit`].
#' Objects produced by the functions [PLNnetwork()], [PLNPCA()], [PLNmixture()] and [PLNLDA()] also enjoy the methods of [PLNfit()] by inheritance.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()], [predict()], [vcov()] and [standard_error()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list-like structure for controlling the fit, see ['PLN_param()'].
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param Theta matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
#' @inherit PLN details
#'
#' @rdname PLNfit
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNfit <- R6Class(
  classname = "PLNfit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## PRIVATE INTERNAL FIELDS
    formula    = NA    , # the formula call for the model as specified by the user
    Theta      = NA    , # regression parameters of the latent layer
    Sigma      = NA    , # covariance matrix of the latent layer
    Omega      = NA    , # precision matrix of the latent layer. Inverse of Sigma
    S          = NA    , # variational parameters for the variances
    M          = NA    , # variational parameters for the means
    Z          = NA    , # matrix of latent variable
    A          = NA    , # matrix of expected counts (under variational approximation)
    Ji         = NA    , # element-wise approximated loglikelihood
    R2         = NA    , # approximated goodness of fit criterion
    optimizer  = list(), # list of links to the functions doing the optimization
    monitoring = list(), # list with optimization monitoring quantities

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torch_elbo = function(data, params) {
      S2 <- torch_multiply(params$S, params$S)
      Z <- data$O + params$M + torch_matmul(data$X, params$Theta)
      res <- .5 * sum(data$w) * torch_logdet(private$torch_Sigma(data, params)) -
        sum(torch_matmul(data$w , data$Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
      res
    },

    torch_Sigma = function(data, params) {
      S2 <- torch_multiply(params$S, params$S)
      Mw <- torch_matmul(torch_diag(torch_sqrt(data$w)), params$M)
      Sigma <- (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(data$w, S2))) / sum(data$w)
      Sigma
    },

    torch_Omega = function(data, params) {
      torch::torch_inverse(params$Sigma)
    },

    torch_vloglik = function(data, params) {
      S2    <- torch_multiply(params$S, params$S)
      Ji <- .5 * self$p - rowSums(.logfactorial(as.matrix(data$Y))) + as.numeric(
        .5 * torch_logdet(params$Omega) +
          torch_sum(data$Y * params$Z - params$A + .5 * torch_log(S2), dim = 2) -
          .5 * torch_sum(torch_matmul(params$M, params$Omega) * params$M + S2 * torch_diag(params$Omega), dim = 2)
      )
      attr(Ji, "weights") <- as.numeric(data$w)
      Ji
    },

    #' @import torch
    torch_optimize = function(Y, X, O, w, init_parameters, configuration) {

      ## Initialization and conversion to torch tensors
      data <- list(
        X = torch_tensor(X),
        Y = torch_tensor(Y),
        O = torch_tensor(O),
        w = torch_tensor(w)
      )

      params <- list(
        Theta = torch_tensor(t(init_parameters$Theta), requires_grad = TRUE),
        M     = torch_tensor(init_parameters$M       , requires_grad = TRUE),
        S     = torch_tensor(init_parameters$S       , requires_grad = TRUE)
      )

      optimizer <- optim_rprop(params, lr = configuration$learning_rate)
      Theta_old <- as.numeric(optimizer$param_groups[[1]]$params$Theta)

      ## Optimization loop
      status <- 5
      objective <- double(length = configuration$maxeval + 1)
      for (iterate in seq.int(configuration$maxeval)) {

        ## Optimization
        optimizer$zero_grad() # reinitialize gradients
        loss <- private$torch_elbo(data, params) # compute current ELBO
        loss$backward()                   # backward propagation
        optimizer$step()                  # optimization

        ## assess convergence
        objective[iterate + 1] <- loss$item()
        Theta_new <- as.numeric(optimizer$param_groups[[1]]$params$Theta)
        delta_f   <- abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1])
        delta_x   <- sum(abs(Theta_old - Theta_new))/sum(abs(Theta_new))
        Theta_old <- Theta_new

        ## display progress
        if (configuration$trace >  1 && (iterate %% 50 == 0))
          cat('\niteration: ', iterate, 'objective', objective[iterate + 1],
              'delta_f'  , round(delta_f, 6), 'delta_x', round(delta_x, 6))

        ## Check for convergence
        if (delta_f < configuration$ftol_rel) status <- 3
        if (delta_x < configuration$xtol_rel) status <- 4
        if (status %in% c(3,4)) {
          objective <- objective[1:iterate + 1]
          break
        }
      }

      params$Sigma <- private$torch_Sigma(data, params)
      params$Omega <- private$torch_Omega(data, params)
      params$Z     <- data$O + params$M + torch_matmul(data$X, params$Theta)
      params$A     <- torch_exp(params$Z + torch_pow(params$S, 2)/2)

      out <- list(
        Theta      = t(as.matrix(params$Theta)),
        Sigma      = as.matrix(params$Sigma),
        Omega      = as.matrix(params$Omega),
        M          = as.matrix(params$M),
        S          = as.matrix(params$S),
        Z          = as.matrix(params$Z),
        A          = as.matrix(params$A),
        Ji         = private$torch_vloglik(data, params),
        monitoring = list(
          objective  = objective,
          iterations = iterate,
          status     = status,
          backend = "torch"
        )
      )
      out
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE METHODS FOR VARIANCE OF THE ESTIMATORS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vcov_sandwich = function(Y, X, Sigma) {
      getMat_iCnTheta <- function(i) {
        a_i   <- as.numeric(private$A[i, ])
        s2_i  <- as.numeric(private$S[i, ]**2)
        # omega <- as.numeric(1/diag(private$Sigma))
        # diag_mat_i <- diag(1/a_i + s2_i^2 / (1 + s2_i * (a_i + omega)))
        diag_mat_i <- diag(1/a_i + .5 * s2_i^2)
        solve(Sigma + diag_mat_i)
      }
      YmA <- Y - private$A
      Dn <- matrix(0, self$d*self$p, self$d*self$p)
      Cn <- matrix(0, self$d*self$p, self$d*self$p)
      for (i in 1:self$n) {
        xxt_i <- tcrossprod(X[i, ])
        Cn <- Cn - kronecker(getMat_iCnTheta(i) , xxt_i) / (self$n)
        Dn <- Dn + kronecker(tcrossprod(YmA[i,]), xxt_i) / (self$n)
      }
      Cn_inv <- solve(Cn)
      (Cn_inv %*% Dn %*% Cn_inv) / (self$n)
    },

    vcov_wald = function(X) {
      fisher <- bdiag(lapply(1:self$p, function(i) {
        crossprod(X, private$A[, i] * X) # t(X) %*% diag(A[, i]) %*% X
      }))
      res <- tryCatch(self$n*Matrix::solve(fisher), error = function(e) {e})
      if (is(res, "error")) {
        warning(paste("Inversion of the Fisher information matrix failed with following error message:",
                      res$message, "Returning NA", sep = "\n"))
        res <- matrix(NA, nrow = self$p, ncol = self$d)
      }
      res
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF PRIVATE METHODS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## CONSTRUCTOR
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description Initialize a [`PLNfit`] model
    #' @importFrom stats lm.wfit lm.fit poisson residuals coefficients runif
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      ## problem dimensions
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)
      ## set up various quantities
      private$formula <- formula # user formula call
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep
      ## initialize the variational parameters
      if (isPLNfit(control$inception)) {
        if (control$trace > 1) cat("\n User defined inceptive PLN model")
        stopifnot(isTRUE(all.equal(dim(control$inception$model_par$Theta), c(p,d))))
        stopifnot(isTRUE(all.equal(dim(control$inception$var_par$M)      , c(n,p))))
        private$Sigma <- control$inception$model_par$Sigma
        private$Theta <- control$inception$model_par$Theta
        private$M     <- control$inception$var_par$M
        private$S     <- sqrt(control$inception$var_par$S2)
      } else {
        if (control$trace > 1) cat("\n Use LM after log transformation to define the inceptive model")
        GLMs <- lapply(1:p, function(j) lm.wfit(covariates, log(1 + responses[,j]), weights, offset = log(1 + exp(offsets[,j]))))
        private$Theta <- do.call(rbind, lapply(GLMs, coefficients))
        private$M     <- do.call(cbind, lapply(GLMs, residuals))
        private$S     <- matrix(1,n,p)
      }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## SETTER METHOD
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description
    #' Update a [`PLNfit`] object
    #' @param M     matrix of variational parameters for the mean
    #' @param S     matrix of variational parameters for the variance
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2    approximate R^2 goodness-of-fit criterion
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param monitoring a list with optimization monitoring quantities
    #' @return Update the current [`PLNfit`] object
    update = function(Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA, Ji=NA, R2=NA, Z=NA, A=NA, monitoring=NA) {
      if (!anyNA(Theta))      private$Theta  <- Theta
      if (!anyNA(Sigma))      private$Sigma  <- Sigma
      if (!anyNA(Omega))      private$Omega  <- Omega
      if (!anyNA(M))          private$M      <- M
      if (!anyNA(S))          private$S      <- S
      if (!anyNA(Z))          private$Z      <- Z
      if (!anyNA(A))          private$A      <- A
      if (!anyNA(Ji))         private$Ji     <- Ji
      if (!anyNA(R2))         private$R2     <- R2
      if (!anyNA(monitoring)) private$monitoring <- monitoring
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## GENERIC OPTIMIZER
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(Y = responses,
                   X = covariates,
                   O = offsets,
                   w = weights,
                   init_parameters = list(Theta = private$Theta, M = private$M, S = private$S),
                   configuration = config)
      optim_out <- do.call(private$optimizer$main, args)
      do.call(self$update, optim_out)
    },

    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, Theta). Intended to position new data in the latent space.
    #' @param Theta Optional fixed value of the regression parameters
    #' @param Sigma variance-covariance matrix of the latent variables
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S2` of variational variances
    #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    optimize_vestep = function(covariates, offsets, responses, weights,
                      Theta = self$model_par$Theta,
                      Omega = self$model_par$Omega,
                      control = PLN_param(backend = "nlopt")) {

      n <- nrow(responses); p <- ncol(responses)
      args <- list(Y = responses,
                   X = covariates,
                   O = offsets,
                   w = weights,
                   ## Initialize the variational parameters with the new dimension of the data
                   init_parameters = list(M = matrix(0, n, p), S = matrix(1, n, p)),
                   Theta = Theta,
                   Omega = Omega,
                   configuration = control$config_optim)
      optim_out <- do.call(private$optimizer$vestep, args)
      optim_out
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Post treatment functions --------------
    #' @description Update R2 field after optimization
    set_R2 = function(responses, covariates, offsets, weights, nullModel = NULL) {
      if (is.null(nullModel)) nullModel <- nullModelPoisson(responses, covariates, offsets, weights)
      loglik <- logLikPoisson(responses, self$latent, weights)
      lmin   <- logLikPoisson(responses, nullModel, weights)
      lmax   <- logLikPoisson(responses, fullModelPoisson(responses, weights), weights)
      private$R2 <- (loglik - lmin) / (lmax - lmin)
    },

    #' @description Experimental: Compute the estimated variance of the coefficient Theta
    #' the true matrix Sigmamust be profided for sandwich estimation at the moment
    #' @param type approximation scheme used, either `wald` (default, variational), `sandwich` (based on MLE theory) or `none`.
    #' @return a sparse matrix with sensible dimension names
    get_vcov_hat = function(type, responses, covariates, Sigma = self$model_par$Sigma) {
      ## compute and store the estimated covariance of the estimator of the parameter Theta
      vcov_hat <-
        switch(type,
               "wald"     = private$vcov_wald(X = covariates),
               "sandwich" = private$vcov_sandwich(Y = responses, X = covariates, Sigma = Sigma),
               "none"     = NULL)

      ## set proper names, use sensible defaults if some names are missing
      rownames(vcov_hat) <- expand.grid(covariates = colnames(covariates),
                                   responses  = colnames(responses)) %>% rev() %>%
        ## Hack to make sure that species is first and varies slowest
        apply(1, paste0, collapse = "_")
      attr(vcov_hat, "name") <- type
      attr(private$Theta, "vcov") <- vcov_hat
    },

    #' @description Update R2, fisher and std_err fields after optimization
    #' @param type approximation scheme used, either `wald` (default, variational), `sandwich` (based on MLE theory) or `none`.
    postTreatment = function(responses, covariates, offsets, weights = rep(1, nrow(responses)), type = c("wald", "sandwich", "none"), nullModel = NULL) {
      type <- match.arg(type)
      ## compute R2
      self$set_R2(responses, covariates, offsets, weights, nullModel)
      ## Set the name of the matrices according to those of the data matrices,
      ## if names are missing, set sensible defaults
      if (is.null(colnames(responses))) colnames(responses) <- paste0("Y", 1:self$p)
      rownames(private$Theta) <- colnames(responses)
      colnames(private$Theta) <- colnames(covariates)
      rownames(private$Sigma) <- colnames(private$Sigma) <- colnames(responses)
      rownames(private$Omega) <- colnames(private$Omega) <- colnames(responses)
      rownames(private$M) <- rownames(private$S) <- rownames(responses)
      colnames(private$S) <- 1:self$q
      if (type != 'none') {
        ## compute and store matrix of standard errors
        self$get_vcov_hat(type, responses, covariates)
      }
    },

    #' @description Predict position, scores or observations of new data.
    #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
    #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    #' @param envir Environment in which the prediction is evaluated
    #' @return A matrix with predictions scores or counts.
    predict = function(newdata, type = c("link", "response"), envir = parent.frame()) {
      type <- match.arg(type)

      ## Extract the model matrices from the new data set with initial formula
      X <- model.matrix(formula(private$formula)[-2], newdata, xlev = attr(private$formula, "xlevels"))
      O <- model.offset(model.frame(formula(private$formula)[-2], newdata))

      ## mean latent positions in the parameter space
      EZ <- tcrossprod(X, private$Theta)
      if (!is.null(O)) EZ <- EZ + O
      EZ <- sweep(EZ, 2, .5 * diag(self$model_par$Sigma), "+")
      colnames(EZ) <- colnames(private$Sigma)

      results <- switch(type, link = EZ, response = exp(EZ))

      attr(results, "type") <- type
      results
    },

    #' @description Predict position, scores or observations of new data, conditionally on the observation of a (set of) variables
    #' @param cond_responses a data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function)
    #' @param newdata a data frame containing the covariates of the sites where to predict
    #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    #' @param var_par Boolean. Should new estimations of the variational parameters of mean and variance be sent back, as attributes of the matrix of predictions. Default to \code{FALSE}.
    #' @param envir Environment in which the prediction is evaluated
    #' @return A matrix with predictions scores or counts.
    predict_cond = function(newdata, cond_responses, type = c("link", "response"), var_par = FALSE, envir = parent.frame()){
      type <- match.arg(type)

      # Checks
      Yc <- as.matrix(cond_responses)
      sp_names <- rownames(self$model_par$Theta)
      if (! any(colnames(cond_responses) %in% sp_names))
        stop("Yc must be a subset of the species in responses")
      if (! nrow(Yc) == nrow(newdata))
        stop("The number of rows of Yc must match the number of rows in newdata")

      # Dimensions and subsets
      n_new <- nrow(Yc)
      cond <- sp_names %in% colnames(Yc)

      ## Extract the model matrices from the new data set with initial formula
      X <- model.matrix(formula(private$formula)[-2], newdata, xlev = attr(private$formula, "xlevels"))
      O <- model.offset(model.frame(formula(private$formula)[-2], newdata))
      if (is.null(O)) O <- matrix(0, n_new, self$p)

      # Compute parameters of the law
      vcov11 <- private$Sigma[cond ,  cond, drop = FALSE]
      vcov22 <- private$Sigma[!cond, !cond, drop = FALSE]
      vcov12 <- private$Sigma[cond , !cond, drop = FALSE]
      prec11 <- solve(vcov11)
      A <- crossprod(vcov12, prec11)
      Sigma21 <- vcov22 - A %*% vcov12

      # Call to VEstep to obtain M1, S1
      VE <- self$optimize_vestep(
              covariates = X,
              offsets    = O[, cond, drop = FALSE],
              responses  = Yc,
              weights    = rep(1, n_new),
              Theta      = self$model_par$Theta[cond, , drop = FALSE],
              Omega      = prec11
          )

      M <- tcrossprod(VE$M, A)
      # S <- map(1:n_new, ~crossprod(sqrt(VE$S[., ]) * t(A)) + Sigma21) %>%
      #   simplify2array()
      S <- map(1:n_new, ~crossprod(VE$S[., ] * t(A)) + Sigma21) %>% simplify2array()

      ## mean latent positions in the parameter space
      EZ <- tcrossprod(X, private$Theta[!cond, , drop = FALSE]) + M
      EZ <- EZ + O[, !cond, drop = FALSE]
      colnames(EZ) <- setdiff(sp_names, colnames(Yc))

      # ! We should only add the .5*diag(S2) term only if we want the type="response"
      if (type == "response") {
        if (ncol(EZ) == 1) {
          EZ <- EZ + .5 * S
        } else {
          EZ <- EZ + .5 * t(apply(S, 3, diag))
        }
      }
      results <- switch(type, link = EZ, response = exp(EZ))
      attr(results, "type") <- type
      if (var_par) {
        attr(results, "M") <- M
        attr(results, "S") <- S
      }
      results
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print functions -----------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    show = function(model = paste("A multivariate Poisson Lognormal fit with", self$vcov_model, "covariance model.\n")) {
      cat(model)
      cat("==================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("==================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $latent, $latent_pos, $var_par, $optim_par\n")
      cat("    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), vcov(), fitted()\n")
      cat("    predict(), predict_cond(), standard_error()\n")
    },

    #' @description User friendly print method
    print = function() { self$show() }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Other functions ----------------
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() {nrow(private$M)},
    #' @field q number of dimensions of the latent space
    q = function() {ncol(private$M)},
    #' @field p number of species
    p = function() {nrow(private$Theta)},
    #' @field d number of covariates
    d = function() {ncol(private$Theta)},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + self$p * (self$p + 1)/2)},
    #' @field model_par a list with the matrices of the model parameters: Theta (covariates), Sigma (covariance), Omega (precision matrix), plus some others depending on the variant)
    model_par  = function() {list(Theta = private$Theta, Sigma = private$Sigma, Omega = private$Omega)},
    #' @field var_par a list with the matrices of the variational parameters: M (means) and S2 (variances)
    var_par    = function() {list(M = private$M, S2 = private$S**2, S = private$S)},
    #' @field optim_par a list with parameters useful for monitoring the optimization
    optim_par  = function() {c(private$monitoring, backend = private$backend)},
    #' @field latent a matrix: values of the latent vector (Z in the model)
    latent     = function() {private$Z},
    #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
    latent_pos = function() {private$M},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$A},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"full"},
    #' @field weights observational weights
    weights     = function() {attr(private$Ji, "weights")},
    #' @field loglik (weighted) variational lower bound of the loglikelihood
    loglik     = function() {sum(self$weights[self$weights > .Machine$double.eps] * private$Ji[self$weights > .Machine$double.eps]) },
    #' @field loglik_vec element-wise variational lower bound of the loglikelihood
    loglik_vec = function() {private$Ji},
    #' @field BIC variational lower bound of the BIC
    BIC        = function() {self$loglik - .5 * log(self$n) * self$nb_param},
    #' @field entropy Entropy of the variational distribution
    entropy    = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(self$var_par$S2)))},
    #' @field ICL variational lower bound of the ICL
    ICL        = function() {self$BIC - self$entropy},
    #' @field R_squared approximated goodness-of-fit criterion
    R_squared  = function() {private$R2},
    #' @field criteria a vector with loglik, BIC, ICL and number of parameters
    criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL)},
    #' @field vcov_coef Approximation of the Variance-Covariance of Theta (experimental)
    vcov_coef  = function() {attr(private$Theta, "vcov")},
    #' @field std_err Approximation of the variance-covariance matrix of model parameters estimates (experimental)
    std_err    = function() {
      if (self$d > 0) {
        stderr <- diag(attr(private$Theta, "vcov")) %>% sqrt %>% matrix(nrow = self$d) %>% t()
        dimnames(stderr) <- dimnames(self$model_par$Theta)
      } else {
        stderr <- NULL
      }
      stderr
    }
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfit_diagonal
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a standard, general framework, with diagonal residual covariance
#'
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
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNfit_diagonal <- R6Class(
  classname = "PLNfit_diagonal",
  inherit = PLNfit,
  public  = list(
    #' @description Initialize a [`PLNfit`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize_diagonal, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_diagonal
    }
  ),
  private = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    torch_elbo = function(data, params) {
      S2 <- torch_pow(params$S, 2)
      Z <- data$O + params$M + torch_matmul(data$X, params$Theta)
      res <- .5 * sum(data$w) * sum(torch_log(private$torch_sigma_diag(data, params))) -
        sum(torch_matmul(data$w , data$Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
      res
    },

    torch_sigma_diag = function(data, params) {
      torch_matmul(data$w, torch_pow(params$M, 2) + torch_pow(params$S,2)) / sum(data$w)
    },

    torch_Sigma = function(data, params) {
      torch_diag(private$torch_sigma_diag(data, params))
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
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + self$p)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"diagonal"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit_diagonal
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfit_spherical
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a standard, general framework, with spherical residual covariance
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#'
#' @rdname PLNfit_spherical
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNfit_spherical <- R6Class(
  classname = "PLNfit_spherical",
  inherit = PLNfit,
  public  = list(
    #' @description Initialize a [`PLNfit`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize_spherical, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_diagonal
    }
  ),
  private = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torch_elbo = function(data, params) {
      S2 <- torch_pow(params$S, 2)
      Z <- data$O + params$M + torch_matmul(data$X, params$Theta)
      res <- .5 * sum(data$w) * self$p * torch_log(private$torch_sigma2(data, params)) -
        sum(torch_matmul(data$w , data$Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
      res
    },

    torch_sigma2 = function(data, params) {
      sum(torch_matmul(data$w, torch_pow(params$M, 2) + torch_pow(params$S, 2))) / (sum(data$w) * self$p)
    },

    torch_Sigma = function(data, params) {
      torch_eye(self$p) * private$torch_sigma2(data, params)
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
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + 1)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"spherical"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit_spherical
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfit_fixedcov
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a standard, general framework, with fixed (inverse) residual covariance
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param config part of the \code{control} argument which configures the optimizer
#'
#' @rdname PLNfit_fixedcov
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNfit_fixedcov <- R6Class(
  classname = "PLNfit_fixedcov",
  inherit = PLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`PLNfit`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$main <- ifelse(control$backend == "nlopt", nlopt_optimize_fixed, private$torch_optimize)
      ## ve step is the same as in the fullly parameterized covariance
      private$Omega <- control$Omega
    },
    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(Y = responses,
                   X = covariates,
                   O = offsets,
                   w = weights,
                   init_parameters = list(Theta = private$Theta, M = private$M, S = private$S, Omega = private$Omega),
                   configuration = config)
      optim_out <- do.call(private$optimizer$main, args)
      do.call(self$update, optim_out)
      private$Sigma <- solve(optim_out$Omega)
    }
  ),
  private = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torch_elbo = function(data, params) {
      S2 <- torch_pow(params$S, 2)
      Z <- data$O + params$M + torch_matmul(data$X, params$Theta)
      res <- sum(data$w) * torch_trace(torch_matmul(private$torch_Sigma(data, params), private$torch_Omega(data, params))) -
        sum(torch_matmul(data$w , data$Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
      res
    },

    torch_Omega = function(data, params) {
      params$Omega <- torch_tensor(private$Omega)
    }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## END OF TORCH METHODS FOR OPTIMIZATION
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # update_loglik = function(weights, Y) {
    #   KY  <- .5 * self$p - rowSums(.logfactorial(Y))
    #   S2 <- private$S**2
    #   Ji <- as.numeric(
    #     .5 * determinant(private$Omega, logarithm = TRUE)$modulus + KY +
    #       rowSums(Y * private$Z - private$A + .5 * log(private$S^2) -
    #                 .5 * ( (private$M %*% private$Omega) * private$M + sweep(private$S^2, 2, diag(private$Omega), '*')))
    #   )
    #   attr(Ji, "weights") <- weights
    #   private$Ji <- Ji
    # },
    #
    # get_objective = function() {
    #   S2 <- S * S
    #   Z <- O + M + torch_matmul(X, Theta)
    #   log_det_Sigma <- switch(configuration$covariance,
    #                           "spherical" = p * log(sum(torch_matmul(w, M * M + S2)) / (w_bar * p)),
    #                           "diagonal"  = sum(torch_log(torch_matmul(w, M * M + S2) / w_bar)),
    #                           { # default value
    #                             Mw <- torch_matmul(torch_diag(torch_sqrt(w)), M)
    #                             Sigma <- (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(w, S2))) / w_bar
    #                             torch_logdet(Sigma)
    #                           })
    #   neg_ELBO <- .5 * w_bar * log_det_Sigma -
    #     sum(torch_matmul(w , Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
    #   neg_ELBO
    # }
  ),
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"fixed"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit_fixedcov
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfit_genetprior
# ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# #' An R6 Class to represent a PLNfit in a standard, general framework, with residual covariance modelling
# #' motivatived by population genetics
# #'
# #' @inherit PLNfit
# #' @rdname PLNfit_genetprior
# #' @importFrom R6 R6Class
# #'
# #' @examples
# #' \dontrun{
# #' data(trichoptera)
# #' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
# #' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
# #' class(myPLN)
# #' print(myPLN)
# #' }
# PLNfit_genetprior <- R6Class(
#   classname = "PLNfit_genetprior",
#   inherit = PLNfit,
#   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   ## PUBLIC MEMBERS ----
#   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   public  = list(
# #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
# optimize = function(responses, covariates, offsets, weights, control) {
#   args <- list(Y = responses,
#                X = covariates,
#                O = offsets,
#                w = weights,
#                init_parameters = list(Theta = private$Theta, M = private$M, S = private$S))
#
#   if (self$vcov_model == "genetic") {
#     args$init_parameters$rho = 0.25
#     args$C <- control$corr_matrix
#   }
#   if (self$vcov_model == "fixed") {
#     args$Omega <- private$Omega
#   }
#
#   if (control$backend == "nlopt")
#     optim_out <- do.call(nlopt_optimizexxx, c(args, list(configuration = control$options_nlopt)))
#   else {
#     ## initialize torch with nlopt
#     optim_out <- self$optimize_nlopt(c(args, list(configuration = control$options_nlopt)))
#     args$init_parameters = list(Theta = optim_out$Theta, M = optim_out$M, S = optim_out$S)
#     optim_out <- self$optimize_torch(c(args, list(configuration = control$options_torch)))
#   }
#
#   private$Theta <- optim_out$Theta
#   private$M     <- optim_out$M
#   private$S     <- optim_out$S
#   private$Z     <- optim_out$Z
#   private$A     <- optim_out$A
#   private$monitoring <- list(iterations = optim_out$iterations, message = status_to_message(optim_out$status))
#   self$update_Sigma(args$w)
#   self$update_loglik(args$w, args$Y)
# },
#     update_Sigma = function(weights) {
#       w_bar <- sum(weights)
#       private$Sigma <- switch(self$vcov_model,
#                               "spherical" = Matrix::Diagonal(self$p, sum(crossprod(weights, private$M^2 + private$S^2)) / (self$p * w_bar)),
#                               "diagonal"  = Matrix::Diagonal(self$p, crossprod(weights, private$M^2 + private$S^2)/ w_bar),
#                               "full"      = (crossprod(private$M, weights * private$M) + diag(as.numeric(crossprod(weights, private$S^2)))) / w_bar,
#                               "fixed"     = solve(private$Omega)
#       )
#       private$Omega <- switch(self$vcov_model,
#                               "fixed"     = private$Omega, solve(private$Sigma)
#                               # "genetic    = private$Omega, solve(private$Sigma)
#       )
#
#       # if (self$vcov_model == "genetic")
#       #   private$psi <- list(sigma2 = optim_out$sigma2, rho = optim_out$rho)
#
#     },
#
#     update_loglik = function(weights, Y) {
#       KY  <- .5 * self$p - rowSums(.logfactorial(Y))
#       S2 <- private$S**2
#       Ji <- as.numeric(
#         .5 * determinant(private$Omega, logarithm = TRUE)$modulus + KY +
#           rowSums(Y * private$Z - private$A + .5 * log(private$S^2) -
#                     .5 * ( (private$M %*% private$Omega) * private$M + sweep(private$S^2, 2, diag(private$Omega), '*')))
#       )
#       attr(Ji, "weights") <- weights
#       private$Ji <- Ji
#     },
#   ),
# active = list(
#   #' @field nb_param number of parameters in the current PLN model
#   nb_param   = function() {as.integer(self$p * self$d + 2)},
#   #' @field vcov_model character: the model used for the residual covariance
#   vcov_model = function() {"genetic"},
# #' @field gen_par a list with two parameters, sigma2 and rho, only used with the genetic covariance model
# gen_par    = function() {private$psi},
# ),
# private = list(
#   psi        = NA, # parameters for genetic model of covariance
# )
#   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   ##  END OF THE CLASS PLNfit_genetprior
#   ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# )
