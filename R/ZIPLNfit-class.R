#' An R6 Class to represent a ZIPLNfit
#'
#' @description The function [ZIPLN()] fit a model which is an instance of a object with class [`ZIPLNfit`].
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()], [predict()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call and used for predictions.
#'
#' @inherit ZIPLN details
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(scRNA)
#' # data subsample: only 100 random cell and the 50 most varying transcript
#' subset <- sample.int(nrow(scRNA), 100)
#' myPLN  <- ZIPLN(counts[, 1:50] ~ 1 + offset(log(total_counts)), subset = subset, data = scRNA)
#' }
ZIPLNfit <- R6Class(
  classname = "ZIPLNfit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Creation functions ----------------
    #' @description
    #' Update a [`ZIPLNfit`] object
    #' @param B     matrix of regression parameters in the Poisson lognormal component
    #' @param B0    matrix of regression parameters in the zero inflated component
    #' @param Pi    Zero inflated probability parameter (either scalar, row-vector, col-vector of matrix)
    #' @param Omega precision matrix of the latent variables
    #' @param Sigma covariance matrix of the latent variables
    #' @param M     matrix of mean vectors for the variational approximation
    #' @param S     matrix of variance parameters for the variational approximation
    #' @param R     matrix of probabilities for the variational approximation
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param monitoring a list with optimization monitoring quantities
    #' @return Update the current [`ZIPLNfit`] object
    update = function(B=NA, B0=NA, Pi=NA, Omega=NA, Sigma=NA, M=NA, S=NA, R=NA, Ji=NA, Z=NA, A=NA, monitoring=NA) {
      if (!anyNA(B))      private$B      <- B
      if (!anyNA(B0))     private$B0     <- B0
      if (!anyNA(Pi))     private$Pi     <- Pi
      if (!anyNA(Omega))  private$Omega  <- Omega
      if (!anyNA(Sigma))  private$Sigma  <- Sigma
      if (!anyNA(M))      private$M      <- M
      if (!anyNA(S))      private$S      <- S
      if (!anyNA(R))      private$R      <- R
      if (!anyNA(Z))      private$Z      <- Z
      if (!anyNA(A))      private$A      <- A
      if (!anyNA(Ji))     private$Ji     <- Ji
      if (!anyNA(monitoring)) private$monitoring <- monitoring
    },

    #' @description Initialize a [`ZIPLNfit`] model
    #' @importFrom stats glm.fit residuals poisson fitted coef
    #' @importFrom pscl zeroinfl
    initialize = function(responses, covariates, offsets, formula, xlevels, control) {
      ## problem dimensions
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)

      ## save the formula call as specified by the user
      private$formula <- formula
      private$xlevels <- xlevels
      private$X       <- covariates
      ## initialize the covariance model
      private$covariance <- control$covariance
      private$ziparam <- control$ziparam
      
      R  <- matrix(0, n, p)
      M  <- matrix(0, n, p)
      B  <- matrix(0, d, p)
      B0 <- matrix(0, d, p)
      Pi <- matrix(0, n, p)
      for (j in 1:p) {
        y = responses[, j]
        if (min(y == 0)) {
          zip_out <- switch(control$ziparam,
            "single" = pscl::zeroinfl(y ~ 0 + covariates | 1  , offset = offsets[, j]),  
            "row"    = pscl::zeroinfl(y ~ 0 + covariates | 1:n, offset = offsets[, j]),  
            "col"    = pscl::zeroinfl(y ~ 0 + covariates | 1  , offset = offsets[, j]),  
            "covar"  = pscl::zeroinfl(y ~ 0 + covariates      , offset = offsets[, j])) # offset only for the count model
          B0[,j] <- coef(zip_out, "zero")
          B[,j]  <- coef(zip_out, "count")
          R[, j] <- predict(zip_out, type = "zero")
          M[,j]  <- residuals(zip_out, type = "pearson") + covariates %*% coef(zip_out, "count")
        } else {
          p_out <- glm(y ~ 0 + covariates, family = 'poisson', offset = offsets[, j])
          B0[,j] <- rep(-10, d)
          B[,j]  <- coef(p_out)
          R[, j] <- 0
          M[,j]  <- residuals(p_out) + covariates %*% coef(p_out)
        }
      }
      
      ## Initialization of the ZI component
      private$R  <- R
      private$Pi <- switch(control$ziparam, 
        "single" = matrix(    mean(R), n, p)              ,
        "row"    = matrix(rowMeans(R), n, p)              ,
        "col"    = matrix(colMeans(R), n, p, byrow = TRUE),
        "covar"  = R)
      private$B0 <- B0
      delta <- 1 * (responses == 0)
      private$zeros <- delta
      
      ## Initialization of the PLN component
      private$B <- B
      private$M <- M
      private$S <- matrix(.1, n, p)
    },

    #' @description Call to the Cpp optimizer and update of the relevant fields
    #' @param control a list for controlling the optimization. See details.
    optimize = function(responses, covariates, offsets, control) {

      args <- list(Y = responses, X = covariates, O = offsets, configuration = control)
      args$init_parameters <- 
        list(Omega = NA, B0 = private$B0, B = private$B, Pi = private$Pi,
             M = private$M, S = private$S, R = private$R)
      optim_out <- do.call(optimize_zi, args)
      
      self$update(
        B      = optim_out$parameters$B,
        B0     = optim_out$parameters$B0,
        Pi     = optim_out$parameters$Pi,
        Omega  = optim_out$parameters$Omega,
        M      = optim_out$parameters$M,
        S      = optim_out$parameters$S,
        R      = optim_out$parameters$R,
        Z      = offsets + optim_out$parameters$M,
        A      = exp(offsets + optim_out$parameters$M + .5 * optim_out$parameters$S^2),
        Ji     = optim_out$vloglik,
        monitoring = list(
          iterations = optim_out$nb_iter,
          message    = optim_out$stop_reason,
          objective  = optim_out$criterion)
      )

      if (is.null(colnames(responses))) colnames(responses) <- paste0("Y", 1:self$p)
      colnames(private$B0) <- colnames(private$B) <- colnames(responses)
      rownames(private$B0) <- rownames(private$B) <- colnames(covariates)
      rownames(private$Omega)  <- colnames(private$Omega) <- colnames(private$Pi) <- colnames(responses)
      rownames(private$M) <- rownames(private$S) <- rownames(private$R) <- rownames(private$Pi) <- rownames(responses)

    },

    # #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, B). Intended to position new data in the latent space.
    # #' @return A list with three components:
    # #'  * the matrix `M` of variational means,
    # #'  * the matrix `S2` of variational variances
    # #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    # VEstep = function(covariates, offsets, responses, control = list()) {
    #
    #   # problem dimension
    #   n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)
    #
    #   ## define default control parameters for optim and overwrite by user defined parameters
    #   control$covariance <- self$vcov_model
    #   control <- ZIPLN_param(control, n, p)
    #
    #   VEstep_optimizer  <-
    #     switch(control$covariance,
    #            "spherical" = cpp_optimize_vestep_spherical,
    #            "diagonal"  = cpp_optimize_vestep_diagonal,
    #            "full"      = cpp_optimize_vestep_full,
    #            "genetic"   = cpp_optimize_vestep_full
    #     )
    #
    #   ## Initialize the variational parameters with the appropriate new dimension of the data
    #   optim_out <- VEstep_optimizer(
    #     list(M = matrix(0, n, p), S = matrix(sqrt(0.1), n, p)),
    #     responses,
    #     covariates,
    #     offsets,
    #     weights,
    #     B = self$model_par$B,
    #     ## Robust inversion using Matrix::solve instead of solve.default
    #     Omega = as(Matrix::solve(Matrix::Matrix(self$model_par$Sigma)), 'matrix'),
    #     control
    #   )
    #
    #   Ji <- optim_out$loglik
    #   attr(Ji, "weights") <- weights
    #
    #   ## output
    #   list(M       = optim_out$M,
    #        S2      = (optim_out$S)**2,
    #        log.lik = setNames(Ji, rownames(responses)))
    # },

    # #' @description Predict position, scores or observations of new data.
    # #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
    # #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    # #' @param envir Environment in which the prediction is evaluated
    # #' @return A matrix with predictions scores or counts.
    # predict = function(newdata, type = c("link", "response"), envir = parent.frame()) {
    #   type = match.arg(type)
    #
    #   ## Extract the model matrices from the new data set with initial formula
    #   X <- model.matrix(formula(private$formula)[-2], newdata, xlev = private$xlevels)
    #   O <- model.offset(model.frame(formula(private$formula)[-2], newdata))
    #
    #   ## mean latent positions in the parameter space
    #   EZ <- X %*% private$B
    #   if (!is.null(O)) EZ <- EZ + O
    #   EZ <- sweep(EZ, 2, .5 * diag(self$model_par$Sigma), "+")
    #   colnames(EZ) <- colnames(private$Sigma)
    #
    #   results <- switch(type, link = EZ, response = exp(EZ))
    #
    #   attr(results, "type") <- type
    #   results
    # },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print functions -----------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    show = function(model = paste("A multivariate Zero Inflated Poisson Lognormal fit with", private$covariance, "covariance model.\n")) {
      cat(model)
      cat("================================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("================================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $var_par, $optim_par, $nb_param\n")
      cat("    $loglik, $BIC, $ICL, $loglik_vec, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), fitted()\n")
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
    formula    = NA, # the formula call for the model as specified by the user
    xlevels    = NA, # factor levels present in the original data, useful for predict() methods.
    X          = NA, # design matrix
    B          = NA, # the model parameters for the covariable effect (PLN part)
    B0         = NA, # the model parameters for the covariate effects ('0'/Bernoulli part)
    Pi         = NA, # the probability parameters for the '0'/Bernoulli part
    Omega      = NA, # the precision matrix
    S          = NA, # the variational parameters for the variances
    M          = NA, # the variational parameters for the means
    Z          = NA, # the matrix of latent variable
    P          = NA, # the matrix of latent variable without covariates effect
    zeros      = NA, # an indicator of the zeros
    A          = NA, # the matrix of expected counts
    R          = NA, # probabilities for being observed
    Ji         = NA, # element-wise approximated loglikelihood
    covariance = NA, # a string describing the covariance model
    ziparam    = NA, # a string describing the ZI parametrisation
    monitoring = NA  # a list with optimization monitoring quantities
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() {nrow(private$M)},
    #' @field p number of species
    p = function() {ncol(private$M)},
    #' @field d number of covariates
    d = function() {nrow(private$B)},
    #' @field latent a matrix: values of the latent vector (Z in the model)
    latent  = function() {private$Z},
    #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
    latent_pos  = function() {private$M - private$X %*% private$B},
    #' @field model_par a list with the matrices of parameters found in the model (B, Sigma, plus some others depending on the variant)
    model_par  = function() {list(B = private$B, B0 = private$B0, Pi = private$Pi, Omega = private$Omega)},
    #' @field var_par a list with two matrices, M and S2, which are the estimated parameters in the variational approximation
    var_par    = function() {list(M = private$M, S2 = private$S^2, R = private$R)},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$R * private$zeros + (1 - private$R) * private$A},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <-  self$p * self$d +
        switch(private$ziparam,
               "single" = 1,
               "row"    = self$n,
               "col"    =  self$p,
               "covar"  =  self$p * self$d) + 
        switch(private$covariance,
          "full" = self$p * (self$p + 1)/2,
          "diagonal" = self$p,
          "spherical" = 1,
          "sparse"    = (sum(private$Omega != 0) - self$p)/2L)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the covariance (either "spherical", "diagonal" or "full")
    vcov_model  = function() {private$covariance},
    #' @field optim_par a list with parameters useful for monitoring the optimization
    optim_par   = function() {private$monitoring},
    #' @field loglik (weighted) variational lower bound of the loglikelihood
    loglik      = function() {sum(private$Ji) },
    #' @field loglik_vec element-wise variational lower bound of the loglikelihood
    loglik_vec  = function() {private$Ji},
    #' @field BIC variational lower bound of the BIC
    BIC         = function() {self$loglik - .5 * log(self$n) * self$nb_param},
    #' @field entropy Entropy of the variational distribution
    entropy     = function() {self$entropy_ZI + self$entropy_PLN},
    #' @field entropy_ZI Entropy of the variational distribution
    entropy_ZI  = function() {-sum(.xlogx(1-private$R)) - sum(.xlogx(private$R))},
    #' @field entropy_PLN Entropy of the Gaussian variational distribution in the PLN component
    entropy_PLN = function() {.5 * (self$n * self$p * log(2*pi*exp(1)) + sum(log(private$S^2)))},
    #' @field ICL variational lower bound of the ICL
    ICL         = function() {self$BIC - self$entropy},
    #' @field criteria a vector with loglik, BIC, ICL and number of parameters
    criteria    = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL)}
  )

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

)

