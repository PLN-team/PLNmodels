#' An R6 Class to represent a ZIPLNfit
#'
#' @description The function [ZIPLN()] fits a model which is an instance of an object with class [`ZIPLNfit`].
#'
#' This class comes with a set of R6 methods, some of which are useful for the end-user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()], [predict()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#'
#' @inherit ZIPLN details
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
#'
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
    #' @param Pi    Zero inflated probability parameter (either scalar, row-vector, col-vector or matrix)
    #' @param Omega precision matrix of the latent variables
    #' @param Sigma covariance matrix of the latent variables
    #' @param M     matrix of mean vectors for the variational approximation
    #' @param S     matrix of standard deviation parameters for the variational approximation
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
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      ## problem dimensions
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates$PLN); d0 <- ncol(covariates$ZI)

      ## save the formula call as specified by the user
      private$formula <- formula
      private$X       <- covariates$PLN
      private$X0      <- covariates$ZI
      ## initialize the covariance model
      private$covariance <- control$covariance
      private$ziparam <- control$ziparam

      R  <- matrix(0, n, p)
      M  <- matrix(0, n, p)
      B  <- matrix(0, d , p)
      B0 <- matrix(0, d0, p)
      ## Feature-wise univariate (ZI)poisson regression as starting point for ZIPLN
      for (j in 1:p) {
        y = responses[, j]
        if (min(y) == 0) {
          suppressWarnings(
            zip_out <- switch(control$ziparam,
              "row"    = pscl::zeroinfl(y ~ 0 + covariates$PLN | 0 + factor(1:n), offset = offsets[, j]),
              "covar"  = pscl::zeroinfl(y ~ 0 + covariates$PLN | 0 + covariates$ZI  , offset = offsets[, j]),
                         pscl::zeroinfl(y ~ 0 + covariates$PLN |               1, offset = offsets[, j])) # offset only for the count model
          )
          B0[,j] <- coef(zip_out, "zero")
          B[,j]  <- coef(zip_out, "count")
          R[, j] <- predict(zip_out, type = "zero")
          M[,j]  <- residuals(zip_out) + covariates$PLN %*% coef(zip_out, "count")
        } else {
          p_out <- glm(y ~ 0 + covariates$PLN, family = 'poisson', offset = offsets[, j])
          B0[,j] <- rep(-10, d)
          B[,j]  <- coef(p_out)
          R[, j] <- 0
          M[,j]  <- residuals(p_out) + covariates$PLN %*% coef(p_out)
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
      private$zeros <- 1 * (responses == 0)

      ## Initialization of the PLN component
      private$B <- B
      private$M <- M
      private$S <- matrix(.1, n, p)

      ## Link to functions performing the optimization
      private$optimizer$B     <- function(M, X) optim_zipln_B_dense(M, X)
      private$optimizer$zi    <- switch(
        control$ziparam,
        "single" = function(R, ...) list(Pi = matrix(    mean(R), nrow(R), p)              , B0 = matrix(NA, d0, p)),
        "row"    = function(R, ...) list(Pi = matrix(rowMeans(R), nrow(R), p)              , B0 = matrix(NA, d0, p)),
        "col"    = function(R, ...) list(Pi = matrix(colMeans(R), nrow(R), p, byrow = TRUE), B0 = matrix(NA, d0, p)),
        "covar"  = optim_zipln_zipar_covar
      )
      private$optimizer$Omega <- optim_zipln_Omega_full

    },

    #' @description Call to the Cpp optimizer and update of the relevant fields
    #' @param control a list for controlling the optimization. See details.
    optimize = function(responses, covariates, offsets, weights, control) {

      data <- list(Y = responses, X = covariates$PLN, X0 = covariates$ZI, O = offsets)
      parameters <-
        list(Omega = NA, B0 = private$B0, B = private$B, Pi = private$Pi,
             M = private$M, S = private$S, R = private$R)

      # Outer loop
      nb_iter <- 0
      criterion   <- numeric(control$maxit_out)
      convergence <- numeric(control$maxit_out)

      vloglik <- -Inf; objective <- Inf
      repeat {

        # Check maxeval
        if (control$maxit_out >= 0 && nb_iter >= control$maxit_out) {
          stop_reason <- "maximum number of iterations reached"
          criterion   <- criterion[1:nb_iter]
          convergence <- convergence[nb_iter]
          break
        }

        ### M Step
        # PLN part
        new_Omega <- private$optimizer$Omega(
          M = parameters$M, X = data$X, B = parameters$B, S = parameters$S
        )
        new_B <- private$optimizer$B(
          M = parameters$M, X = data$X
        )

        # ZI part
        optim_new_zipar <- private$optimizer$zi(
          R = parameters$R, init_B0 = parameters$B0, X0 = data$X0, config = control
        )
        new_B0 <- optim_new_zipar$B0
        new_Pi <- optim_new_zipar$Pi

        ### VE Step
        # ZI part
        new_R <- optim_zipln_R(
          Y = data$Y, X = data$X, O = data$O, M = parameters$M, S = parameters$S, Pi = new_Pi
        )
        # PLN part
        new_M <- optim_zipln_M(
          init_M = parameters$M,
          Y = data$Y, X = data$X, O = data$O, R = new_R, S = parameters$S, B = new_B, Omega = new_Omega,
          configuration = control
        )$M
        new_S <- optim_zipln_S(
          init_S = parameters$S,
          O = data$O, M = new_M, R = new_R, B = new_B, diag_Omega = diag(new_Omega),
          configuration = control
        )$S

        # Check convergence
        new_parameters <- list(
          Omega = new_Omega, B = new_B, B0 = new_B0, Pi = new_Pi,
          R = new_R, M = new_M, S = new_S
        )
        nb_iter <- nb_iter + 1

        vloglik <- zipln_vloglik(
          data$Y, data$X, data$O, new_Pi, new_Omega, new_B, new_R, new_M, new_S
        )

        criterion[nb_iter] <- new_objective <- -sum(vloglik)
        convergence[nb_iter]  <- abs(new_objective - objective)/abs(new_objective)

        objective_converged <-
          (objective - new_objective) <= control$ftol_out |
          (objective - new_objective)/abs(new_objective) <= control$ftol_out

        parameters_converged <- parameter_list_converged(
          parameters, new_parameters,
          xtol_abs = control$xtol_abs, xtol_rel = control$xtol_rel
        )

        if (parameters_converged | objective_converged) {
          parameters <- new_parameters
          stop_reason <- "converged"
          criterion   <- criterion[1:nb_iter]
          convergence <- convergence[1:nb_iter]
          break
        }

        parameters <- new_parameters
        objective  <- new_objective
      }

      self$update(
        B      = parameters$B,
        B0     = parameters$B0,
        Pi     = parameters$Pi,
        Omega  = parameters$Omega,
        Sigma =  tryCatch(Matrix::solve(symmpart(parameters$Omega)), error = function(e) {e}),
        M      = parameters$M,
        S      = parameters$S,
        R      = parameters$R,
        Z      = offsets + parameters$M,
        A      = exp(offsets + parameters$M + .5 * parameters$S^2),
        Ji     = vloglik,
        monitoring = list(
          iterations  = nb_iter,
          message     = stop_reason,
          objective   = criterion,
          convergence = convergence)
      )

      ### TODO: Should be in post-treatment
      if (is.null(colnames(responses))) colnames(responses) <- paste0("Y", 1:self$p)
      colnames(private$B0) <- colnames(private$B) <- colnames(responses)
      rownames(private$B0) <- rownames(private$B) <- colnames(covariates)
      rownames(private$Omega)  <- colnames(private$Omega) <- colnames(private$Pi) <- colnames(responses)
      dimnames(private$Sigma)  <- dimnames(private$Omega)
      rownames(private$M) <- rownames(private$S) <- rownames(private$R) <- rownames(private$Pi) <- rownames(responses)

    },

    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S, R) and corresponding log likelihood values for fixed model parameters (Sigma, B, B0). Intended to position new data in the latent space.
    #' @param B Optional fixed value of the regression parameters in the PLN component
    #' @param B0 Optional fixed value of the regression parameters in the ZI component
    #' @param Omega inverse variance-covariance matrix of the latent variables
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S` of variational standard deviations
    #'  * the matrix `R` of variational ZI probabilities
    #'  * the vector `Ji` of (variational) log-likelihood of each new observation
    #'  * a list `monitoring` with information about convergence status
    optimize_vestep = function(covariates, offsets, responses, weights,
                               B = self$model_par$B,
                               B0 = self$model_par$B0,
                               Omega = self$model_par$Omega,
                               control = ZIPLN_param(backend = "nlopt")$config_optim) {
      n <- nrow(responses)
      data <- list(Y = responses, X = covariates$PLN, X0 = covariates$ZI, O = offsets)
      parameters <-
        list(M = matrix(0, n, self$p), S = matrix(0.1, n, self$p), R = matrix(0, n, self$p))

      # Outer loop
      nb_iter <- 0
      criterion   <- numeric(control$maxit_out)
      convergence <- numeric(control$maxit_out)

      vloglik <- -Inf; objective <- Inf

      repeat {

        # Check maxeval
        if (control$maxit_out >= 0 && nb_iter >= control$maxit_out) {
          stop_reason <- "maximum number of iterations reached"
          criterion   <- criterion[1:nb_iter]
          convergence <- convergence[nb_iter]
          break
        }

        Pi <- private$optimizer$zi(
          R = parameters$R, init_B0 = B0, X0 = data$X0, config = config_default_nlopt
        )$Pi

        # VE Step
        new_R <- optim_zipln_R(
          Y = data$Y, X = data$X, O = data$O, M = parameters$M, S = parameters$S, Pi = Pi
        )
        new_M <- optim_zipln_M(
          init_M = parameters$M,
          Y = data$Y, X = data$X, O = data$O, R = new_R, S = parameters$S, B = B, Omega = Omega,
          configuration = control
        )$M
        new_S <- optim_zipln_S(
          init_S = parameters$S,
          O = data$O, M = new_M, R = new_R, B = B, diag_Omega = diag(Omega),
          configuration = control
        )$S
        # Check convergence
        new_parameters <- list(R = new_R, M = new_M, S = new_S)
        nb_iter <- nb_iter + 1

        vloglik <- zipln_vloglik(
          data$Y, data$X, data$O, Pi, Omega, B, new_R, new_M, new_S
        )

        criterion[nb_iter] <- new_objective <- -sum(vloglik)
        convergence[nb_iter]  <- abs(new_objective - objective)/abs(new_objective)

        objective_converged <-
          (objective - new_objective) <= control$ftol_out |
          (objective - new_objective)/abs(new_objective) <= control$ftol_out

        parameters_converged <- parameter_list_converged(
          parameters, new_parameters,
          xtol_abs = control$xtol_abs, xtol_rel = control$xtol_rel
        )

        ## Update parameters
        parameters <- new_parameters
        objective <- new_objective

        ## End outer loop in case of convergence
        if (parameters_converged | objective_converged) {
          stop_reason <- "converged"
          criterion   <- criterion[1:nb_iter]
          convergence <- convergence[1:nb_iter]
          break
        }

      }

      list(
        M      = parameters$M,
        S      = parameters$S,
        R      = parameters$R,
        Ji     = vloglik,
        monitoring = list(
          iterations  = nb_iter,
          message     = stop_reason,
          objective   = criterion,
          convergence = convergence)
      )
    },

    #' @description Predict position, scores or observations of new data.
    #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are returned.
    #' @param responses Optional data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function), assuming the interest is in testing the model.
    #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    #' @param level Optional integer value the level to be used in obtaining the predictions. Level zero corresponds to the population predictions (default if `responses` is not provided) while level one (default) corresponds to predictions after evaluating the variational parameters for the new data.
    #' @param envir Environment in which the prediction is evaluated
    #'
    #' @details
    #' Note that `level = 1` can only be used if responses are provided,
    #' as the variational parameters can't be estimated otherwise. In the absence of responses, `level` is ignored and the fitted values are returned
    #' @return A matrix with predictions scores or counts.
    predict = function(newdata, responses = NULL, type = c("link", "response"), level = 1, envir = parent.frame()) {

      ## Ignore everything if newdata is not provided
      if (missing(newdata)) {
        return(self$fitted)
      }

      n_new <- nrow(newdata)
      ## Set level to 0 (to bypass VE step) if responses are not provided
      if (is.null(responses)) {
        level <- 0
      }

      terms <- .extract_terms_zi(as.formula(private$formula))

      ## Extract the model matrices from the new data set with initial formula
      # PLN part
      X <- model.matrix(terms$PLN[-2], newdata, xlev = attr(private$formula, "xlevels")$PLN)
      # ZI part
      if (!is.null(terms$ZI)) {
        X0 <- model.matrix(terms$ZI, newdata, xlev = attr(private$formula, "xlevels")$ZI)
      } else {
        X0 <- matrix(NA,0,0)
      }

      O <- model.offset(model.frame(terms$PLN[-2], newdata))
      if (is.null(O)) O <- matrix(0, n_new, self$p)

      ## Optimize M and S if responses are provided,
      if (level == 1) {
        VE <- self$optimize_vestep(
          covariates = list(PLN = X, ZI = X0),
          offsets    = O,
          responses  = as.matrix(responses),
          weights    = rep(1, n_new),
          B          = private$B,
          B0         = private$B0,
          Omega      = private$Omega
        )
        R <- VE$R
        M <- VE$M
        S2 <- VE$S^2
      } else {
        # otherwise set R to Pi, M to XB and S2 to diag(Sigma)
        R <- private$Pi[1:nrow(newdata), ]
        M <- X %*% private$B
        S2 <- matrix(diag(private$Sigma), nrow = n_new, ncol = self$p, byrow = TRUE)
      }

      ## mean latent positions in the parameter space (covariates/offset only)
      EZ <- O + M
      rownames(EZ) <- rownames(newdata)
      colnames(EZ) <- colnames(private$Sigma)

      type <- match.arg(type)
      results <- switch(
        type,
        link = EZ,
        response = R + (1 - R) * exp(EZ + .5 * S2),

      )
      attr(results, "type") <- type
      results
    },

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
    X          = NA, # design matrix for the PLN component
    X0         = NA, # design matrix for the ZI component
    B          = NA, # the model parameters for the covariate effects (PLN part)
    B0         = NA, # the model parameters for the covariate effects (ZI part)
    Pi         = NA, # the probability parameters for the ZI part
    Omega      = NA, # the precision matrix
    Sigma      = NA, # the covariance matrix
    S          = NA, # the variational parameters for the standard deviations
    M          = NA, # the variational parameters for the means
    Z          = NA, # the matrix of latent variable
    P          = NA, # the matrix of latent variable without covariates effect
    zeros      = NA, # an indicator of the zeros
    A          = NA, # the matrix of expected counts
    R          = NA, # probabilities for being observed
    Ji         = NA, # element-wise approximated loglikelihood
    covariance = NA, # a string describing the covariance model
    ziparam    = NA, # a string describing the ZI model (single, col, row, covar)
    optimizer  = list(), # list of links to the functions doing the optimization
    monitoring = list()  # list with optimization monitoring quantities
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples/sites
    n = function() {nrow(private$M)},
    #' @field q number of dimensions of the latent space
    q = function() {ncol(private$M)},
    #' @field p number of variables/species
    p = function() {ncol(private$B)},
    #' @field d number of covariates in the PLN part
    d = function() {nrow(private$B)},
    #' @field d0 number of covariates in the ZI part
    d0 = function() {nrow(private$B0)},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      as.integer(
        self$p * self$d + self$p * (self$p + 1L)/2L +
          switch(private$ziparam,
                 "single" = 1L,
                 "row"    = self$n,
                 "col"    =  self$p,
                 "covar"  =  self$p * self$d)
      )
    },
    #' @field model_par a list with the matrices of parameters found in the model (B, Sigma, plus some others depending on the variant)
    model_par  = function() {list(B = private$B, B0 = private$B0, Pi = private$Pi, Omega = private$Omega, Sigma = private$Sigma)},
    #' @field var_par a list with two matrices, M and S2, which are the estimated parameters in the variational approximation
    var_par    = function() {list(M = private$M, S2 = private$S^2, S = private$S, R = private$R)},
    #' @field optim_par a list with parameters useful for monitoring the optimization
    optim_par   = function() {private$monitoring},
    #' @field latent a matrix: values of the latent vector (Z in the model)
    latent  = function() {private$Z},
    #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
    latent_pos  = function() {private$M - private$X %*% private$B},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$R * private$zeros + (1 - private$R) * private$A},
    #' @field vcov_model character: the model used for the covariance (either "spherical", "diagonal", "full" or "sparse")
    vcov_model  = function() {private$covariance},
    #' @field zi_model character: the model used for the zero inflation (either "single", "row", "col" or "covar")
    zi_model  = function() {private$ziparam},
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
  ##  END OF THE CLASS ZIPLNfit
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZIPLNfit_diagonal ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a ZIPLNfit in a standard, general framework, with diagonal residual covariance
#'
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
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "diagonal"))
#' class(myPLN)
#' print(myPLN)
#' }
ZIPLNfit_diagonal <- R6Class(
  classname = "ZIPLNfit_diagonal",
  inherit = ZIPLNfit,
  public  = list(
    #' @description Initialize a [`ZIPLNfit_diagonal`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$Omega <- optim_zipln_Omega_diagonal
    }
  ),
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <-  self$p * self$d + self$p +
        switch(private$ziparam,
               "single" = 1L,
               "row"    = self$n,
               "col"    =  self$p,
               "covar"  =  self$p * self$d)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"diagonal"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ZIPLNfit_diagonal
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZIPLNfit_spherical  ##########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a ZIPLNfit in a standard, general framework, with spherical residual covariance
#'
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
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "spherical"))
#' class(myPLN)
#' print(myPLN)
#' }
ZIPLNfit_spherical <- R6Class(
  classname = "ZIPLNfit_spherical",
  inherit = ZIPLNfit,
  public  = list(
    #' @description Initialize a [`ZIPLNfit_spherical`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$Omega <- optim_zipln_Omega_spherical
    }
  ),
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <-  self$p * self$d + 1L +
        switch(private$ziparam,
               "single" = 1L,
               "row"    = self$n,
               "col"    =  self$p,
               "covar"  =  self$p * self$d)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"spherical"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit_spherical
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZIPLNfit_fixed  #############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a ZIPLNfit in a standard, general framework, with fixed (inverse) residual covariance
#'
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
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera,
#'     control = ZIPLN_param(Omega = diag(ncol(trichoptera$Abundance))))
#' class(myPLN)
#' print(myPLN)
#' }
ZIPLNfit_fixed <- R6Class(
  classname = "ZIPLNfit_fixed",
  inherit = ZIPLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`ZIPLNfit_fixed`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$Omega <- control$Omega
      private$optimizer$Omega <- function(M, X, B, S) {private$Omega}
    }
  ),
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <-  self$p * self$d +
        switch(private$ziparam,
               "single" = 1L,
               "row"    = self$n,
               "col"    =  self$p,
               "covar"  =  self$p * self$d)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"fixed"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ZIPLNfit_fixed
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZIPLNfit_sparse  #############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a ZIPLNfit in a standard, general framework, with sparse inverse residual covariance
#'
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
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control=  ZIPLN_param(penalty = 0.2))
#' class(myPLN)
#' print(myPLN)
#' }
ZIPLNfit_sparse <- R6Class(
  classname = "ZIPLNfit_sparse",
  inherit = ZIPLNfit,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`ZIPLNfit_fixed`] model
    #' @importFrom glassoFast glassoFast
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$Omega <-
        function(M, X, B, S) {
          glassoFast( crossprod(M - X %*% B)/self$n + diag(colMeans(S * S), self$p, self$p), rho = control$penalty )$wi
        }
    }
  ),
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {
      res <-  self$p * self$d + (sum(private$Omega != 0) - self$p)/2L +
        switch(private$ziparam,
               "single" = 1L,
               "row"    = self$n,
               "col"    =  self$p,
               "covar"  =  self$p * self$d)
      as.integer(res)
    },
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"sparse"}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ZIPLNfit_sparse
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
