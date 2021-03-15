#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNmixture}} produces a collection of models which are instances of object with class \code{PLNmixturefit}.
#' A \code{PLNmixturefit} (say, with k components) is itself a collection of k \code{PLNfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for ...
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param control a list for controlling the optimization. See details.
#' @param clusters the dimensions of the successively fitted models
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param xlevels named listed of factor levels included in the models, extracted from the formula in the upper-level call #'
#' @param cluster the number of clusters of the current model
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#'
#' @include PLNfit-class.R
#'
#' @importFrom R6 R6Class
#' @importFrom purrr map2_dbl map2
#' @importFrom purrr map_int map_dbl
#' @seealso The function \code{\link{PLNmixture}}, the class \code{\link[=PLNmixturefamily]{PLNmixturefamily}}
PLNmixturefit <-
  R6Class(classname = "PLNmixturefit",
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      formula    = NA, # the formula call for the model as specified by the user
      xlevels    = NA, # factor levels present in the original data, useful for predict() methods.
      covariance = NA, # a string describing the covariance model
      comp       = NA, # list of mixture components (PLNfit)
      tau        = NA, # posterior probabilities of cluster belonging
      monitoring = NA, # a list with optimization monitoring quantities
      Theta      = NA, # the model parameters for the covariable
      mix_up     = function(var_name) {
        private$comp %>%
          map(function(C) C$weights * eval(str2expression(paste0('C$', var_name)))) %>%
          reduce(`+`)
      },
      #' @description Optimize a the
      optimize_covariates = function(Y, X, O){

        M  <- private$comp %>%  map("var_par") %>% map("M")
        S2 <- private$comp %>%  map("var_par") %>% map("S2")
        if(private$covariance == "spherical") S2 <- map(S2, ~outer(as.numeric(.x), rep(1, self$p) ))

        mu <- private$comp %>%  map(coef) %>% map(~outer(rep(1, self$n), as.numeric(.x)))

        Ak_tilde <- list(M, S2, mu) %>%
          purrr::pmap(function(M_k, S2_k, mu_k) exp(O + mu_k + M_k + .5 * S2_k))

        Tk <- asplit(private$tau, 2)

        objective_and_gradient <- function(theta) {

          Theta <- matrix(theta, ncol(X), self$p)

          XB <- X %*% Theta
          Ak <- map(Ak_tilde, ~.x * exp(XB))

          list("objective" = sum(purrr::map2_dbl(Tk, Ak, ~sum (t(.x) %*% (.y - Y * XB) ))),
               "gradient"  = Reduce('+', purrr::map2(Tk, Ak, ~ as.numeric(t(diag(.x) %*% X) %*% (.y - Y)))))
        }

        opts <- list("algorithm"="NLOPT_LD_CCSAQ", "xtol_rel"=1.0e-4)
        out_optim <- nloptr::nloptr(c(private$Theta), objective_and_gradient,
                                    opts = opts)
        private$Theta <- matrix(out_optim$solution, ncol(X), self$p)
        invisible(out_optim)
      }
    ),
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public  = list(
      #' @description Initialize a [`PLNmixturefit`] model
      #'@param posteriorProb matrix ofposterior probability for cluster belonging
      initialize = function(responses, covariates, offsets, posteriorProb, formula, xlevels, control) {
        private$tau   <- posteriorProb
        private$comp  <- vector('list', ncol(posteriorProb))
        private$Theta <- matrix(0, ncol(covariates), ncol(responses))
        private$formula <- formula
        private$xlevels <- xlevels
        private$covariance <- control$covariance

        ## Initializing the mixture components (only intercept of group mean)
        mu_k  <- setNames(matrix(1, self$n, ncol = 1), 'Intercept')
        for (k_ in seq.int(ncol(posteriorProb)))
          private$comp[[k_]] <- PLNfit$new(responses, mu_k, offsets, posteriorProb[, k_], NULL, NULL, control)

      },
      #' @description Optimize a [`PLNmixturefit`] model
      optimize = function(responses, covariates, offsets, control) {

        ## The intercept term will serve as the mean in each group/component
        intercept <- matrix(1, nrow(responses), ncol = 1)
        ## We make a copy of the offset, for accounting for fixed
        ## covariates effects during the alternative algorithm
        offsets_ <- offsets

        ## ===========================================
        ## INITIALISATION
        cond <- FALSE; iter <- 1
        objective   <- numeric(control$maxit_out); objective[iter]   <- Inf
        convergence <- numeric(control$maxit_out); convergence[iter] <- NA
        ## ===========================================
        ## OPTIMISATION
        while (!cond) {
          iter <- iter + 1
          if (control$trace > 1) cat("", iter)

          ## ---------------------------------------------------
          ## M - STEP
          ## UPDATE Theta, THE MATRIX OF REGRESSION COEFFICIENTS
          if (ncol(covariates) > 0) {
            private$optimize_covariates(responses, covariates, offsets_)
            offsets <- offsets_ + covariates %*% private$Theta
          }
          ## UPDATE THE MIXTURE MODEL VIA OPTIMIZATION OF PLNmixture
          for (k in seq.int(self$k))
            self$components[[k]]$optimize(responses, intercept, offsets, private$tau[, k], control)

          ## ---------------------------------------------------
          ## E - STEP
          ## UPDATE THE POSTERIOR PROBABILITIES
          if (self$k > 1) { # only needed when at least 2 components!
            private$tau <-
              sapply(private$comp, function(comp) comp$loglik_vec) %>% # Jik
              sweep(2, log(self$mixtureParam), "+") %>% # computation in log space
              apply(1, .softmax) %>%        # exponentiation + normalization with soft-max
              t() %>% .check_boundaries()   # bound away probabilities from 0/1
          }

          ## Assess convergence
          objective[iter]   <- -self$loglik
          convergence[iter] <- abs(objective[iter-1] - objective[iter]) /abs(objective[iter])
          if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE

        }

        ## ===========================================
        ## OUTPUT
        ## formatting parameters for output
        private$monitoring = list(objective        = objective[2:iter],
                                  convergence      = convergence[2:iter],
                                  outer_iterations = iter-1)
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Prediction methods --------------------
      #' @description Predict group of new samples
      #' @param newdata A data frame in which to look for variables, offsets and counts with which to predict.
      #' @param type The type of prediction required. The default `posterior` are posterior probabilities for each group ,
      #'  `response` is the group with maximal posterior probability and `latent` is the averaged latent coordinate (without
      #'  offset and nor covariate effects),
      #'  with weights equal to the posterior probabilities.
      #' @param prior User-specified prior group probabilities in the new data. The default uses a uniform prior.
      #' @param control a list for controlling the optimization. See [PLN()] for details.
      #' @param envir Environment in which the prediction is evaluated
      predict = function(newdata,
                         type = c("posterior", "response", "position"),
                         prior = matrix(rep(1/self$k, self$k), nrow(newdata), self$k, byrow = TRUE),
                         control = list(), envir = parent.frame()) {

        type  <- match.arg(type)

        ## Extract the model matrices from the new data set with initial formula
        args <- extract_model(call("PLNmixture", formula = private$formula, data = newdata, xlev = private$xlevels), envir)
        n_new <- nrow(args$Y)

        ## Sanity checks
        stopifnot(ncol(args$X) == self$d)
        stopifnot(ncol(args$O) == self$p, ncol(args$Y) == self$p)
        stopifnot(nrow(prior)  == nrow(args$X), ncol(prior) == self$k)
        stopifnot(nrow(args$X) == n_new)

        ## The intercept term will serve as the mean in each group/component
        intercept <- matrix(1, n_new, ncol = 1)
        ## Initialize the posteriorProbabilities
        tau <- prior
        ## We make a copy of the offset, for accounting for fixed
        ## covariates effects during the alternative algorithm
        if (ncol(args$X) > 1) {
          args$O <- args$O + args$X %*% private$Theta
        }

        ## ===========================================
        ## INITIALISATION
        cond <- FALSE; iter <- 1
        control <- PLNmixture_param(control, n_new, self$p)
        objective   <- numeric(control$maxit_out); objective[iter]   <- Inf
        convergence <- numeric(control$maxit_out); convergence[iter] <- NA

        ## ===========================================
        ## OPTIMISATION
        while(!cond) {
          iter <- iter + 1
          if (control$trace > 1) cat("", iter)

          ## ---------------------------------------------------
          ## VE step of each component
          ve_step <- list(self$k)
          for (k in seq.int(self$k)) {
            ve_step[[k]] <- self$components[[k]]$VEstep(intercept, args$O, args$Y, tau[, k], control = control)
          }

          ## E - STEP
          ## UPDATE THE POSTERIOR PROBABILITIES
          if (self$k > 1) { # only needed when at least 2 components!
            tau <-
              sapply(ve_step, function(comp) comp$log.lik) %>% # Jik
              sweep(2, log(colMeans(tau)), "+") %>% # computation in log space
              apply(1, .softmax) %>%        # exponentiation + normalization with soft-max
              t() %>% .check_boundaries()   # bound away probabilities from 0/1
          }

          ## Assess convergence
          J_ik <- sapply(ve_step, function(comp) comp$log.lik)
          J_ik[tau <= .Machine$double.eps] <- 0
          rowSums(tau * J_ik) - rowSums(.xlogx(tau)) + tau %*% log(colMeans(tau))
          objective[iter]   <- -sum(J_ik)
          convergence[iter] <- abs(objective[iter-1] - objective[iter]) /abs(objective[iter])
          if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE

        }

        switch(type,
          "posterior" = {rownames(tau) <- rownames(newdata); tau},
          "response"  = as.factor(setNames(apply(tau, 1, which.max), rownames(newdata))),
          "position"  = {
            latent_pos <- array(0, dim = c(nrow(args$X), self$k, self$p))
            for (k in seq.int(self$k)) {
              latent_pos[ , k, ] <- ve_step[[k]]$M + rep(1, nrow(args$X)) %o% self$group_means[, k]
            }
            res <- apply(latent_pos * tau %o% rep(1, self$p), c(1, 3), sum)
            rownames(res) <- rownames(newdata)
            res
            }
          )
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------
      #' @description Plot the matrix of expected mean counts (without offsets, without covariate effects) reordered according the inferred clustering
      #' @param plot logical. Should the plot be displayed or sent back as [`ggplot`] object
      #' @param main character. A title for the plot.  An hopefully appropriate title will be used by default.
      #' @param log_scale logical. Should the color scale values be log-transform before plotting? Default is \code{TRUE}.
      #' @return a [`ggplot`] graphic
      plot_clustering_data = function(main = "Expected counts reorder by clustering", plot = TRUE, log_scale = TRUE) {
        M  <- private$mix_up('var_par$M')
        S2 <- switch(private$covariance, "spherical" = private$mix_up('var_par$S2') %*% rbind(rep(1, ncol(M))), private$mix_up('var_par$S2'))
        mu <- self$posteriorProb %*% t(self$group_means)
        A  <- exp(mu + M + .5 * S2)
        p <- plot_matrix(A, 'samples', 'variables', self$memberships, log_scale)
        if (plot) print(p)
        invisible(p)
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------
      #' @description Plot the individual map of a PCA performed on the latent coordinates, where individuals are colored according to the memberships
      #' @param plot logical. Should the plot be displayed or sent back as [`ggplot`] object
      #' @param main character. A title for the plot. An hopefully appropriate title will be used by default.
      #' @return a [`ggplot`] graphic
      plot_clustering_pca = function(main = "Clustering labels in Individual Factor Map", plot = TRUE) {
        mu <- self$posteriorProb %*% t(self$group_means) + private$mix_up('var_par$M')
        svdM <- svd(scale(mu, TRUE, FALSE), nv = 2)
        .scores <- data.frame(t(t(svdM$u[, 1:2]) * svdM$d[1:2]))
        colnames(.scores) <- paste("a",1:2,sep = "")
        .scores$labels <- as.factor(self$memberships)
        .scores$names <- rownames(self$components[[1]]$var_par$M)
        eigen.val <- svdM$d^2
        .percent_var <- round(eigen.val[1:2]/sum(eigen.val),4)
        axes_label <- paste(paste("axis",1:2), paste0("(", round(100* .percent_var,3)[1:2], "%)"))
        p <- get_ggplot_ind_map(.scores, axes_label, main)
        if (plot) print(p)
        invisible(p)
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Post treatment --------------------
      #' @description Update fields after optimization
      #' @param weights an optional vector of observation weights to be used in the fitting process.
      postTreatment = function(responses, covariates, offsets, weights, nullModel) {

        ## restoring the full design matrix (group means + covariates)
        mu_k <- setNames(matrix(1, self$n, ncol = 1), 'Intercept')
        offsets <- offsets + covariates %*% private$Theta

        for (k_ in seq.int(ncol(private$tau)))
          self$components[[k_]]$postTreatment(
            responses,
            mu_k,
            offsets,
            private$tau[,k_],
            nullModel = nullModel,
            type = "none"
          )
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Print methods ---------------------
      #' @description User friendly print method
      show = function() {
        cat("Poisson Lognormal mixture model with",self$k,"components and", self$vcov_model,"covariances.\n")
        cat("* Useful fields\n")
        cat("    $posteriorProb, $memberships, $mixtureParam, $group_means\n")
        cat("    $model_par, $latent, $latent_pos, $optim_par\n")
        cat("    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria\n")
        cat("    $component[[i]] (a PLNfit with associated methods and fields)\n")
        cat("* Useful S3 methods\n")
        cat("    print(), coef(), sigma(), fitted(), predict() \n")
      },
      #' @description User friendly print method
      print = function() self$show()
    ),
    active = list(
      #' @field n number of samples
      n = function() {nrow(private$tau)},
      #' @field p number of dimensions of the latent space
      p = function() {ncol(self$model_par$Theta)},
      #' @field k number of components
      k = function() {length(private$comp)},
      #' @field d number of covariates
      d = function() {nrow(private$Theta)},
      #' @field components components of the mixture (PLNfits)
      components    = function(value) {if (missing(value)) return(private$comp) else private$comp <- value},
      #' @field latent a matrix: values of the latent vector (Z in the model)
      latent = function() {private$mix_up('latent')},
      #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
      latent_pos = function() {private$mix_up('var_par$M') + self$posteriorProb %*% t(self$group_means)},
      #' @field posteriorProb matrix ofposterior probability for cluster belonging
      posteriorProb = function(value) {if (missing(value)) return(private$tau) else private$tau <- value},
      #' @field memberships vector for cluster index
      memberships = function(value) {apply(private$tau, 1, which.max)},
      #' @field mixtureParam vector of cluster proportions
      mixtureParam  = function() {colMeans(private$tau)},
      #' @field optim_par a list with parameters useful for monitoring the optimization
      optim_par  = function() {private$monitoring},
      #' @field nb_param number of parameters in the current PLN model
      nb_param  = function() {(self$k-1) + self$p * self$d + sum(map_int(self$components, 'nb_param'))},
      #' @field entropy_clustering Entropy of the variational distribution of the cluster (multinomial)
      entropy_clustering = function() {-sum(.xlogx(private$tau))},
      #' @field entropy_latent Entropy of the variational distribution of the latent vector (Gaussian)
      entropy_latent = function() {
        .5 * (sum(map_dbl(private$comp, function(component) {
          sum( diag(component$weights) %*% log(component$var_par$S2) * ifelse(component$vcov_model == "spherical", self$p, 1) )
          })) + self$n * self$p * log(2*pi*exp(1)))
      },
      #' @field entropy Full entropy of the variational distribution (latent vector + clustering)
      entropy = function() {self$entropy_latent + self$entropy_clustering},
      #' @field loglik variational lower bound of the loglikelihood
      loglik = function() {sum(self$loglik_vec)},
      #' @field loglik_vec element-wise variational lower bound of the loglikelihood
      loglik_vec = function() {
        J_ik <- sapply(private$comp, function(comp_) comp_$loglik_vec)
        J_ik[private$tau <= .Machine$double.eps] <- 0
        rowSums(private$tau * J_ik) - rowSums(.xlogx(private$tau)) + private$tau %*% log(self$mixtureParam)
        },
      #' @field BIC variational lower bound of the BIC
      BIC       = function() {self$loglik - .5 * log(self$n) * self$nb_param},
      #' @field ICL variational lower bound of the ICL (include entropy of both the clustering and latent distributions)
      ICL       = function() {self$BIC - self$entropy},
      #' @field R_squared approximated goodness-of-fit criterion
      R_squared = function() {sum(self$mixtureParam * map_dbl(self$components, "R_squared"))},
      #' @field criteria a vector with loglik, BIC, ICL, and number of parameters
      criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL)},
      #' @field model_par a list with the matrices of parameters found in the model (Theta, Sigma, Mu and Pi)
      model_par  = function() {list(Theta = private$Theta,
                                    Sigma = private$comp %>% map('model_par') %>% map('Sigma'),
                                    Mu = self$group_means,
                                    Pi = self$mixtureParam)},
      #' @field vcov_model character: the model used for the covariance (either "spherical", "diagonal" or "full")
      vcov_model = function() {private$covariance},
      #' @field fitted a matrix: fitted values of the observations (A in the model)
      fitted = function() {private$mix_up('fitted')},
      #' @field group_means a matrix of group mean vectors in the latent space.
      group_means = function() {
        self$components %>%
          map(function(C) C$model_par$Theta)  %>%
          setNames(paste0("group_", 1:self$k)) %>% as.data.frame()
      }
    )
)
