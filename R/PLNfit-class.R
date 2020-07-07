#' An R6 Class to represent a PLNfit in a standard, general framework
#'
#' @description The function [PLN()] fit a model which is an instance of a object with class [`PLNfit`].
#' Objects produced by the functions [PLNnetwork()], [PLNPCA()] and [PLNLDA()] also enjoy the methods of [PLNfit()] by inheritance.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()],
#' [predict()], [vcov()] and [standard_error()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param model model used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call and used for predictions.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param type approximation scheme to compute the fisher information matrix. Either `wald` (default) or `louis`. \code{type = "louis"} results in smaller confidence intervals.
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
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description
    #' Update a [`PLNfit`] object
    #' @param Theta matrix of regression matrix
    #' @param Sigma variance-covariance matrix of the latent variables
    #' @param M     matrix of mean vectors for the variational approximation
    #' @param S2    matrix of variance vectors for the variational approximation
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2    approximate R^2 goodness-of-fit criterion
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param monitoring a list with optimization monitoring quantities
    #' @return Update the current [`PLNfit`] object
    update = function(Theta=NA, Sigma=NA, M=NA, S2=NA, Ji=NA, R2=NA, Z=NA, A=NA, monitoring=NA) {
      if (!anyNA(Theta))      private$Theta  <- Theta
      if (!anyNA(Sigma))      private$Sigma  <- Sigma
      if (!anyNA(M))          private$M      <- M
      if (!anyNA(S2))         private$S2     <- S2
      if (!anyNA(Z))          private$Z      <- Z
      if (!anyNA(A))          private$A      <- A
      if (!anyNA(Ji))         private$Ji     <- Ji
      if (!anyNA(R2))         private$R2     <- R2
      if (!anyNA(monitoring)) private$monitoring <- monitoring
    },

    #' @description Initialize a [`PLNfit`] model
    #' @importFrom stats lm.wfit lm.fit poisson residuals coefficients runif
    ## TODO: Once "set" is supported by Roxygen go back to external definition using
    ## PLNfit$set("public", "initialize", { ... })
    initialize = function(responses, covariates, offsets, weights, model, xlevels, control) {
      ## problem dimensions
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)

      ## save the formula call as specified by the user
      private$model      <- model
      private$xlevels    <- xlevels
      ## initialize the covariance model
      private$covariance <- control$covariance
      private$optimizer  <-
        switch(control$covariance,
               "spherical" = optim_spherical,
               "diagonal"  = optim_diagonal ,
               "full"      = optim_full     ,
               "rank"      = optim_rank     ,
               "sparse"    = optim_sparse
        )

      if (isPLNfit(control$inception)) {
        if (control$trace > 1) cat("\n User defined inceptive PLN model")
        stopifnot(isTRUE(all.equal(dim(control$inception$model_par$Theta), c(p,d))))
        stopifnot(isTRUE(all.equal(dim(control$inception$var_par$M)      , c(n,p))))
        private$Theta <- control$inception$model_par$Theta
        private$M     <- control$inception$var_par$M
        private$S2    <- control$inception$var_par$S2
        private$Sigma <- control$inception$model_par$Sigma
        private$Ji    <- control$inception$loglik_vec
      } else {
        # if (control$trace > 1) cat("\n Use GLM Poisson to define the inceptive model")
        # LMs   <- lapply(1:p, function(j) glm.fit(covariates, responses[,j], weights, offset =  offsets[,j], family = poisson(), intercept = FALSE))
        if (control$trace > 1) cat("\n Use LM after log transformation to define the inceptive model")
        LMs   <- lapply(1:p, function(j) lm.wfit(covariates, log(1 + responses[,j]), weights, offset =  offsets[,j]) )
        private$Theta <- do.call(rbind, lapply(LMs, coefficients))
        residuals     <- do.call(cbind, lapply(LMs, residuals))
        private$M     <- residuals
        private$S2    <- matrix(0.1, n, ifelse(control$covariance == "spherical", 1, p))
        if (control$covariance == "spherical") {
          private$Sigma <- diag(sum(residuals^2)/(n*p), p, p)
        } else  if (control$covariance == "diagonal") {
          private$Sigma <- diag(diag(crossprod(residuals)/n), p, p)
        } else  {
          private$Sigma <- crossprod(residuals)/n + diag(colMeans(private$S2), nrow = p)
        }
      }

    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Optimizers ----------------------------
    #' @description Call to the C++ optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, control) {
      optim_out <- private$optimizer(
        c(private$Theta, private$M, sqrt(private$S2)),
        responses,
        covariates,
        offsets,
        weights,
        control
      )

      Ji <- optim_out$loglik
      attr(Ji, "weights") <- weights
      self$update(
        Theta      = optim_out$Theta,
        Sigma      = optim_out$Sigma,
        M          = optim_out$M,
        S2         = (optim_out$S)**2,
        Z          = optim_out$Z,
        A          = optim_out$A,
        Ji         = Ji,
        monitoring = list(
          iterations = optim_out$iterations,
          status     = optim_out$status,
          message    = statusToMessage(optim_out$status))
      )
    },

    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, Theta). Intended to position new data in the latent space.
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S` of variational variances
    #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    VEstep = function(covariates, offsets, responses, weights, control = list()) {

      # problem dimension
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)

      ## define default control parameters for optim and overwrite by user defined parameters
      control$covariance <- self$vcov_model
      control <- PLN_param(control, n, p, d)

      VEstep_optimizer  <-
        switch(control$covariance,
               "spherical" = VEstep_PLN_spherical,
               "diagonal"  = VEstep_PLN__diagonal,
               "full"      = VEstep_PLN_full
        )

      optim_out <- VEstep_optimizer(
        c(private$M, sqrt(private$S2)),
        responses,
        covariates,
        offsets,
        weights,
        Theta = self$model_par$Theta,
        ## Robust inversion using Matrix::solve instead of solve.default
        Omega = as(Matrix::solve(Matrix::Matrix(self$model_par$Sigma)), 'matrix'),
        control
      )

      ## output
      list(M       = optim_out$M,
           S2      = (optim_out$S)**2,
           log.lik = setNames(optim_out$loglik, rownames(responses)))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Post treatment functions --------------
    #' @description Update R2 field after optimization
    set_R2 = function(responses, covariates, offsets, weights, nullModel = NULL) {
      if (is.null(nullModel)) nullModel <- nullModelPoisson(responses, covariates, offsets, weights)
      loglik <- logLikPoisson(responses, self$latent_pos(covariates, offsets), weights)
      lmin   <- logLikPoisson(responses, nullModel, weights)
      lmax   <- logLikPoisson(responses, fullModelPoisson(responses, weights), weights)
      private$R2 <- (loglik - lmin) / (lmax - lmin)
    },

    #' @description Safely compute the fisher information matrix (FIM)
    #' @param X design matrix used to compute the FIM
    #' @return a sparse matrix with sensible dimension names
    compute_fisher = function(type = c("wald", "louis"), X = NULL) {
      type <- match.arg(type)
      A <- private$A
      if (type == "louis") {
        ## TODO check how to adapt for PLNPCA
        ## A = A + A \odot A \odot (exp(S) - 1_{n \times p})
        A <- A + A * A * (exp(self$var_par$S) - 1)
      }
      if (anyNA(A)) {
        warning("Something went wrong during model fitting!!\nMatrix A has missing values.")
        result <- bdiag(lapply(1:self$p, function(i) {diag(NA, nrow = self$d)}))
      } else {
        result <- bdiag(lapply(1:self$p, function(i) {
          ## t(X) %*% diag(A[, i]) %*% X
          crossprod(X, A[, i] * X)
        }))
      }
      ## set proper names, use sensible defaults if some names are missing
      element.names <- expand.grid(covariates = colnames(private$Theta),
                                   species    = rownames(private$Theta)) %>% rev() %>%
        ## Hack to make sure that species is first and varies slowest
        apply(1, paste0, collapse = "_")
      rownames(result) <- element.names
      result
    },

    #' @description Compute univariate standard error for coefficients of Theta from the FIM
    #' @return a matrix of standard deviations.
    #' @importFrom Matrix diag solve
    compute_standard_error = function() {
      if (self$d > 0) {
        ## self$fisher$mat : Fisher Information matrix I_n(\Theta) = n * I(\Theta)
        ## safe inversion using Matrix::solve and Matrix::diag and error handling
        out <- tryCatch(Matrix::diag(Matrix::solve(self$fisher$mat)),
                        error = function(e) {e})
        if (is(out, "error")) {
          warning(paste("Inversion of the Fisher information matrix failed with following error message:",
                        out$message,
                        "Returning NA",
                        sep = "\n"))
          stderr <- matrix(NA, nrow = self$p, ncol = self$d)
        } else {
          stderr <- out %>% sqrt %>% matrix(nrow = self$d) %>% t()
        }
        dimnames(stderr) <- dimnames(self$model_par$Theta)
      } else {
        stderr <- NULL
      }
      stderr
    },

    #' @description Update R2, fisher and std_err fields after optimization
    postTreatment = function(responses, covariates, offsets, weights = rep(1, nrow(responses)), type = c("wald", "louis"), nullModel = NULL) {
      ## compute R2
      self$set_R2(responses, covariates, offsets, weights, nullModel)
      ## Set the name of the matrices according to those of the data matrices,
      ## if names are missing, set sensible defaults
      if (is.null(colnames(responses))) colnames(responses) <- paste0("Y", 1:self$p)
      rownames(private$Theta) <- colnames(responses)
      colnames(private$Theta) <- colnames(covariates)
      rownames(private$Sigma) <- colnames(private$Sigma) <- colnames(responses)
      rownames(private$M) <- rownames(private$S2) <- rownames(responses)
      ## compute and store Fisher Information matrix
      type <- match.arg(type)
      private$FIM <- self$compute_fisher(type, X = covariates)
      private$FIM_type <- type
      ## compute and store matrix of standard errors
      private$.std_err <- self$compute_standard_error()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Helper functions ----------------------
    #' @description Compute matrix of latent positions, noted as Z in the model. Used to compute the likelihood or for data visualization
    #' @return a n x q matrix of latent positions.
    latent_pos = function(covariates, offsets) {
      latentPos <- private$M + tcrossprod(covariates, private$Theta) + offsets
      latentPos
    },

    #' @description Predict position, scores or observations of new data.
    #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
    #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    #' @param envir Environment in which the prediction is evaluated
    #' @return A matrix with predictions scores or counts.
    predict = function(newdata, type = c("link", "response"), envir = parent.frame()) {
      type = match.arg(type)

      ## Extract the model matrices from the new data set with initial formula
      X <- model.matrix(formula(private$model)[-2], newdata, xlev = private$xlevels)
      O <- model.offset(model.frame(formula(private$model)[-2], newdata))

      ## mean latent positions in the parameter space
      EZ <- tcrossprod(X, private$Theta)
      if (!is.null(O)) EZ <- EZ + O
      EZ <- sweep(EZ, 2, .5 * diag(self$model_par$Sigma), "+")
      colnames(EZ) <- colnames(private$Sigma)

      results <- switch(type, link = EZ, response = exp(EZ))

      attr(results, "type") <- type
      results
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print functions -----------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    show = function(model = paste("A multivariate Poisson Lognormal fit with", private$covariance, "covariance model.\n")) {
      cat(model)
      cat("==================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("==================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $latent, $var_par, $optim_par\n")
      cat("    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), vcov(), fitted(), predict(), standard_error()\n")
    },

    #' @description User friendly print method
    print = function() { self$show() }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Other functions ----------------
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    model      = NA, # the formula call for the model as specified by the user
    xlevels    = NA, # factor levels present in the original data, useful for predict() methods.
    Theta      = NA, # the model parameters for the covariable
    Sigma      = NA, # the covariance matrix
    S2         = NA, # the variational parameters for the variances
    M          = NA, # the variational parameters for the means
    Z          = NA, # the matrix of latent variable
    A          = NA, # the matrix of expected counts
    R2         = NA, # approximated goodness of fit criterion
    Ji         = NA, # element-wise approximated loglikelihood
    FIM        = NA, # Fisher information matrix of Theta, computed using of two approximation scheme
    FIM_type   = NA, # Either "wald" or "louis". Approximation scheme used to compute FIM
    .std_err   = NA, # element-wise standard error for the elements of Theta computed
    # from the Fisher information matrix
    covariance = NA, # a string describing the covariance model
    optimizer  = NA, # link to the function that performs the optimization
    monitoring = NA  # a list with optimization monitoring quantities
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
    #' @field model_par a list with the matrices of parameters found in the model (Theta, Sigma, plus some others depending on the variant)
    model_par  = function() {list(Theta = private$Theta, Sigma = private$Sigma)},
    #' @field fisher Variational approximation of the Fisher Information matrix
    fisher     = function() {list(mat = private$FIM, type = private$FIM_type) },
    #' @field std_err Variational approximation of the variance-covariance matrix of model parameters estimates.
    std_err    = function() {private$.std_err },
    #' @field var_par a list with two matrices, M and S2, which are the estimated parameters in the variational approximation
    var_par    = function() {list(M = private$M, S2 = private$S2)},
    #' @field latent a matrix: values of the latent vector (Z in the model)
    latent     = function() {private$Z},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$A},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <- self$p * self$d + switch(private$covariance, "full" = self$p * (self$p + 1)/2, "diagonal" = self$p, "spherical" = 1)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the covariance (either "spherical", "diagonal" or "full")
    vcov_model = function() {private$covariance},
    #' @field optim_par a list with parameters useful for monitoring the optimization
    optim_par  = function() {private$monitoring},
    #' @field loglik (weighted) variational lower bound of the loglikelihood
    loglik     = function() {sum(attr(private$Ji, "weights") * private$Ji) },
    #' @field loglik_vec element-wise variational lower bound of the loglikelihood
    loglik_vec = function() {private$Ji},
    #' @field BIC variational lower bound of the BIC
    BIC        = function() {self$loglik - .5 * log(self$n) * self$nb_param},
    #' @field entropy Entropy of the variational distribution
    entropy    = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(private$S2)) * ifelse(private$covariance == "spherical", self$q, 1))},
    #' @field ICL variational lower bound of the ICL
    ICL        = function() {self$BIC - self$entropy},
    #' @field R_squared approximated goodness-of-fit criterion
    R_squared  = function() {private$R2},
    #' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
    criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL, R_squared = self$R_squared)}
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

)


