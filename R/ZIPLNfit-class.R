## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## CLASS ZIPLNfit -----
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a ZIPLNfit
#'
#' @description The function [ZIPLN()] fits a model which is an instance of an object with class [`ZIPLNfit`].
#'
#' This class comes with a set of R6 methods, some of which are useful for the end-user and exported as S3 methods.
#' See the documentation for [coef()], [sigma()], [predict()].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @param data a named list used internally to carry the data matrices
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
    #' @importFrom tidyr replace_na
    initialize = function(data, control) {
      ## problem dimensions
      n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); d0 <- ncol(data$X0)

      ## save the formula call as specified by the user
      private$formula <- data$formula
      private$X       <- data$X
      private$X0      <- data$X0
      ## initialize the covariance model
      private$covariance <- control$covariance
      private$ziparam    <- control$ziparam
      if (isZIPLNfit(control$inception)) {
        private$R  <- control$inception$var_par$R
        private$M  <- control$inception$var_par$M
        private$S  <- control$inception$var_par$S
        private$B  <- control$inception$model_par$B
        private$B0 <- control$inception$model_par$B0
      } else {
        R  <- matrix(0, n, p)
        M  <- matrix(0, n, p)
        B  <- matrix(0, d , p)
        B0 <- matrix(0, d0, p)
        ## Feature-wise univariate (ZI)poisson regression as starting point for ZIPLN
        for (j in 1:p) {
          y = data$Y[, j]
          if (min(y) == 0) {
            suppressWarnings(
              zip_out <- switch(control$ziparam,
                  "row"    = pscl::zeroinfl(y ~ 0 + data$X | 0 + factor(1:n), offset = data$O[, j]),
                  "covar"  = pscl::zeroinfl(y ~ 0 + data$X | 0 + data$X0    , offset = data$O[, j]),
                             pscl::zeroinfl(y ~ 0 + data$X |               1, offset = data$O[, j])) # offset only for the count model
            )
            B0[,j] <- replace_na(coef(zip_out, "zero") , -10)
            B[,j]  <- replace_na(coef(zip_out, "count"), 0)
            R[, j] <- replace_na(predict(zip_out, type = "zero"), sum(y == 0) / n)
            M[,j]  <- replace_na(residuals(zip_out), 0) + data$X %*% coef(zip_out, "count")
          } else {
            p_out  <- glm(y ~ 0 + data$X, family = 'poisson', offset = data$O[, j])
            B0[,j] <- rep(-10, d)
            B[,j]  <- replace_na(coef(p_out), 0)
            R[, j] <- sum(y == 0) / n
            M[,j]  <- replace_na(residuals(p_out), 0) + data$X %*% coef(p_out)
          }
        }

        ## Initialization of the ZI component
        private$R  <- R
        private$B0 <- B0
        ## Initialization of the PLN component
        private$B <- B
        private$M <- M
        private$S <- matrix(.1, n, p)
      }
      private$Pi <- switch(control$ziparam,
                           "single" = matrix(    mean(private$R), n, p)              ,
                           "row"    = matrix(rowMeans(private$R), n, p)              ,
                           "col"    = matrix(colMeans(private$R), n, p, byrow = TRUE),
                           "covar"  = private$R)
      private$zeros <- 1 * (data$Y == 0)

      ## Link to functions performing the optimization
      private$optimizer$B     <- function(M, X) optim_zipln_B_dense(M, X)
      private$optimizer$zi    <- switch(
        control$ziparam,
        "single" = function(R, ...) list(Pi = matrix(    mean(R), nrow(R), p)              , B0 = matrix(NA, d0, p)),
        "row"    = function(R, ...) list(Pi = matrix(rowMeans(R), nrow(R), p)              , B0 = matrix(NA, d0, p)),
        "col"    = function(R, ...) list(Pi = matrix(colMeans(R), nrow(R), p, byrow = TRUE), B0 = matrix(NA, d0, p)),
        "covar"  = optim_zipln_zipar_covar
      )
      private$optimizer$R <- ifelse(control$config_optim$approx_ZI, optim_zipln_R_var, optim_zipln_R_exact)
      private$optimizer$Omega <- optim_zipln_Omega_full

    },

    #' @description Call to the Cpp optimizer and update of the relevant fields
    #' @param control a list for controlling the optimization. See details.
    optimize = function(data, control) {

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
        new_R <- private$optimizer$R(Y = data$Y, X = data$X, O = data$O, M = parameters$M, S = parameters$S, Pi = new_Pi, B = new_B)

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
        Z      = data$O + parameters$M,
        A      = exp(data$O + parameters$M + .5 * parameters$S^2),
        Ji     = vloglik,
        monitoring = list(
          iterations  = nb_iter,
          message     = stop_reason,
          objective   = criterion,
          convergence = convergence)
      )

      ### TODO: Should be in post-treatment
      colnames_Y <- colnames(data$Y)
      if (is.null(colnames_Y)) colnames_Y <- paste0("Y", 1:self$p)
      colnames(private$B0) <- colnames(private$B) <- colnames_Y
      rownames(private$B)  <- colnames(data$X)
      rownames(private$B0) <- colnames(data$X0)
      rownames(private$Omega)  <- colnames(private$Omega) <- colnames(private$Pi) <- colnames_Y
      dimnames(private$Sigma)  <- dimnames(private$Omega)
      rownames(private$M) <- rownames(private$S) <- rownames(private$R) <- rownames(private$Pi) <- rownames(data$Y)

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
    optimize_vestep = function(data,
                               B = self$model_par$B,
                               B0 = self$model_par$B0,
                               Omega = self$model_par$Omega,
                               control = ZIPLN_param(backend = "nlopt")$config_optim) {
      n <- nrow(data$Y)
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
        new_R <- private$optimizer$R(
          Y = data$Y, X = data$X, O = data$O, M = parameters$M, S = parameters$S, Pi = Pi, B = B
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

    #' @description Predict position, scores or observations of new data. See [predict.ZIPLNfit()] for the S3 method and additional details
    #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
    #' @param responses Optional data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function), assuming the interest in in testing the model.
    #' @param type Scale used for the prediction. Either `"link"` (default, predicted positions in the latent space), `"response"` (predicted average counts, accounting for zero-inflation) or `"deflated"` (predicted average counts, not accounting for zero-inflation and using only the PLN part of the model).
    #' @param level Optional integer value the level to be used in obtaining the predictions. Level zero corresponds to the population predictions (default if `responses` is not provided) while level one (default) corresponds to predictions after evaluating the variational parameters for the new data.
    #' @param envir Environment in which the prediction is evaluated
    #'
    #' @return A matrix with predictions scores or counts.
    predict = function(newdata, responses = NULL, type = c("link", "response", "deflated"), level = 1, envir = parent.frame()) {

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
          data = list(
            Y = as.matrix(responses), X = X, X0 = X0, O = O, w = rep(1, n_new)
          ),
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
        link     = EZ,
        response = (1 - R) * exp(EZ + .5 * S2),
        deflated = exp(EZ + .5 * S2),
      )
      attr(results, "type") <- type
      results
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print functions -----------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    show = function(model = paste("A multivariate Zero Inflated Poisson Lognormal fit with", self$vcov_model, "covariance model.\n")) {
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
    #' @field nb_param_zi number of parameters in the ZI part of the model
    nb_param_zi   = function() {
      as.integer(switch(private$ziparam,
                        "single" = 1L,
                        "row"    = self$n,
                        "col"    =  self$p,
                        "covar"  =  self$p * self$d))
    },
    #' @field nb_param_pln number of parameters in the PLN part of the model
    nb_param_pln   = function() {
      as.integer(self$p * self$d + self$p * (self$p + 1L) / 2L)
    },
    #' @field nb_param number of parameters in the ZIPLN model
    nb_param   = function() {
      self$nb_param_zi + self$nb_param_pln
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
    fitted     = function() {(1 - private$R) * private$A},
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
#' @param data a named list used internally to carry the data matrices
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
    initialize = function(data, control) {
      super$initialize(data, control)
      private$optimizer$Omega <- optim_zipln_Omega_diagonal
    }
  ),
  active = list(
    #' @field nb_param_pln number of parameters in the PLN part of the current model
    nb_param_pln   = function() {
      as.integer(self$p * self$d + self$p)
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
#' @param data a named list used internally to carry the data matrices
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
    initialize = function(data, control) {
      super$initialize(data, control)
      private$optimizer$Omega <- optim_zipln_Omega_spherical
    }
  ),
  active = list(
    #' @field nb_param_pln number of parameters in the PLN part of the current model
    nb_param_pln   = function() {
      as.integer(self$p * self$d + 1L)
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
#' @param data a named list used internally to carry the data matrices
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
    initialize = function(data, control) {
      super$initialize(data, control)
      private$Omega <- control$Omega
      private$optimizer$Omega <- function(M, X, B, S) {private$Omega}
    }
  ),
  active = list(
    #' @field nb_param_pln number of parameters in the PLN part of the current model
    nb_param_pln   = function() {
      as.integer(self$p * self$d + 0L)
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
#' @param data a named list used internally to carry the data matrices
#' @param control a list for controlling the optimization. See details.
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' # See other examples in function ZIPLN
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control=  ZIPLN_param(penalty = 1))
#' class(myPLN)
#' print(myPLN)
#' plot(myPLN)
#' }
ZIPLNfit_sparse <- R6Class(
  classname = "ZIPLNfit_sparse",
  inherit = ZIPLNfit,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    lambda = NA, # the sparsity tuning parameter
    rho    = NA  # the p x p penalty weight
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`ZIPLNfit_fixed`] model
    #' @importFrom glassoFast glassoFast
    initialize = function(data, control) {
      super$initialize(data, control)
      ## Default for penalty weights (if not already set)
      if (is.null(control$penalty_weights)) control$penalty_weights <- matrix(1, self$p, self$p)
      stopifnot(isSymmetric(control$penalty_weights), all(control$penalty_weights >= 0))
      if (!control$penalize_diagonal) diag(control$penalty_weights) <- 0
      private$lambda <- control$penalty
      private$rho    <- control$penalty_weights
      private$optimizer$Omega <-
        function(M, X, B, S) {
          glassoFast( crossprod(M - X %*% B)/self$n + diag(colMeans(S * S), self$p, self$p),
                      rho = private$lambda * private$rho )$wi
        }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @importFrom Matrix Matrix
    #' @return a square matrix of size `ZIPLNfit_sparse$n`
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
      .plot_network(self$latent_network(match.arg(type)),
                    type            = match.arg(type),
                    output          = match.arg(output),
                    edge.color      = edge.color,
                    remove.isolated = remove.isolated,
                    node.labels     = node.labels,
                    layout          = layout,
                    plot            = plot)
    }
  ),
  active = list(
    #' @field penalty the global level of sparsity in the current model
    penalty         = function() {private$lambda},
    #' @field penalty_weights a matrix of weights controlling the amount of penalty element-wise.
    penalty_weights = function() {private$rho},
    #' @field n_edges number of edges if the network (non null coefficient of the sparse precision matrix)
    n_edges         = function() {sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0)},
    #' @field nb_param_pln number of parameters in the PLN part of the current model
    nb_param_pln   = function() {
      as.integer(self$p * self$d + self$n_edges + self$p)
    },
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"sparse"},
    #' @field pen_loglik variational lower bound of the l1-penalized loglikelihood
    pen_loglik      = function() {self$loglik - private$lambda * sum(abs(private$Omega))},
    #' @field EBIC variational lower bound of the EBIC
    EBIC      = function() {self$BIC - .5 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$p*(self$p - 1)/self$n_edges), 0)},
    #' @field density proportion of non-null edges in the network
    density   = function() {mean(self$latent_network("support"))},
    #' @field criteria a vector with loglik, penalized loglik, BIC, EBIC, ICL, R_squared, number of parameters, number of edges and graph density
    criteria  = function() {data.frame(super$criteria, n_edges = self$n_edges, EBIC = self$EBIC, pen_loglik = self$pen_loglik, density = self$density)}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS ZIPLNfit_sparse
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)
