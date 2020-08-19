#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNmixture}} produces a collection of models which are instances of object with class \code{PLNmixturefit}.
#' A \code{PLNmixturefit} (say, with k components) is itself a collection of k \code{PLNfit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for ...
#'
#' @param cluster the number of clusters of the current model
#' @param components a list with cluster component element, each of whom is a \code{PLNfit}.
#' @param model_par a list with the matrices associated with the final estimated parameters of the mixture model: Theta (covariates), Sigma (latent covariance), mu (vector of means/centers) and pi (vector of cluster proportions)
#' @param posteriorProbabilities matrix of posterior probabilities of class belonging
#' @param mixtureParam vector of cluster proportions
#' @param loglik variational lower bound of the loglikelihood
#' @param BIC variational lower bound of the BIC
#' @param ICL variational lower bound of the ICL
#' @param R_squared approximated goodness-of-fit criterion
#'
#' @include PLNfit-class.R
#'
#' @importFrom R6 R6Class
#' @importFrom parallel mclapply
#' @importFrom purrr map_int map_dbl
#' @seealso The function \code{\link{PLNmixture}}, the class \code{\link[=PLNmixturefamily]{PLNmixturefamily}}
PLNmixturefit <-
  R6Class(classname = "PLNmixturefit",
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      comp       = NA, # list of mixture components (PLNfit)
      tau        = NA, # posterior probabilities of cluster belonging
      monitoring = NA,  # a list with optimization monitoring quantities
      mix_up     = function(var_name) {
        Reduce("+",
           Map(function(pi, comp) {
             pi * eval(str2expression(paste0('comp$', var_name)))
          }, self$mixtureParam, self$components)
        )
      }
    ),
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public  = list(
      #' @description Initialize a [`PLNmixturefit`] model
      initialize = function(responses, covariates, offsets, posteriorProb, model, xlevels, control) {
        private$tau  <- posteriorProb
        private$comp <- vector('list', ncol(posteriorProb))
        for (k_ in seq.int(ncol(posteriorProb)))
          private$comp[[k_]] <- PLNfit$new(responses, covariates, offsets, posteriorProb[, k_], model, xlevels, control)
      },
      #' @description Optimize a [`PLNmixturefit`] model
      optimize = function(responses, covariates, offsets, control) {
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
            ## UPDATE THE MIXTURE MODEL VIA OPTIMIZATION OF PLNmixture
            # parallel::mclapply(seq.int(self$k), function(k_){
            #     self$components[[k_]]$optimize(responses, covariates, offsets, private$tau[, k_], control)
            # }, mc.cores = control$cores)

            for (k_ in seq.int(self$k))
              self$components[[k_]]$optimize(responses, covariates, offsets, private$tau[, k_], control)

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

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------
      #' @description Plot the matrix of mean counts (without offsets, without covariate effects) reordered according the inferred clustering
      #' @param plot logical. Should the plot be displayed or sent back as [`ggplot`] object
      #' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
      #' @param log_scale logical. Should the color scale values be log-transform before plotting? Default is \code{TRUE}.
      #' @return a [`ggplot`] graphic
      plot_clustering_data = function(main = "Expected counts reorder by clustering", plot = TRUE, log_scale = TRUE) {
        M  <- self$var_par$M
        S2 <- self$var_par$S2
        A  <- exp(M + .5 * S2 %*% rbind(rep(1,ncol(M))))
        p <- plot_matrix(A, 'samples', 'variables', self$memberships, log_scale)
        if (plot) print(p)
        invisible(p)
      },
      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Graphical methods -----------------
      #' @description Plot the individual map of a PCA performed on the latent coordinate, where individuals are colored according to the memberships
      #' @param plot logical. Should the plot be displayed or sent back as [`ggplot`] object
      #' @param main character. A title for the single plot (individual or variable factor map). If NULL (the default), an hopefully appropriate title will be used.
      #' @return a [`ggplot`] graphic
      plot_clustering_pca = function(main = "Clustering labels in Individual Factor Map", plot = TRUE) {
        svdM <- svd(self$var_par$M, nv = 2)
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
      postTreatment = function(responses, covariates, offsets, weights, nullModel) {
        for (k_ in seq.int(ncol(private$tau)))
          self$components[[k_]]$postTreatment(
            responses,
            covariates,
            offsets,
            private$tau[,k_],
            nullModel = nullModel
          )
      },
      show = function() {
        cat("Poisson Lognormal mixture model with",self$k,"components.\n")
        cat("* check fields $posteriorProb, $memberships, $mixtureParam and $components\n")
        cat("* check methods $plot_clustering_data, $plot_clustering_pca\n")
        cat("* each $component[[i]] is a PLNfit with expacted methods and fields\n")
      },
      print = function() self$show()
    ),
    active = list(
      #' @field n number of samples
      n = function() {nrow(private$tau)},
      #' @field k number of components
      k = function() {length(private$comp)},
      #' @field components components of the mixture (PLNfits)
      components    = function(value) {if (missing(value)) {return(private$comp)} else {private$comp <- value}},
      #' @field posteriorProb matrix ofposterior probability for cluster belonging
      posteriorProb = function() {private$tau},
      #' @field meberships vector for cluster index
      memberships   = function(value) {apply(private$tau, 1, which.max)},
      #' @field mixtureParam vector of cluster proporitions
      mixtureParam  = function() {colMeans(private$tau)},
      #' @field optim_par a list with parameters useful for monitoring the optimization
      optim_par  = function() {private$monitoring},
      #' @field nb_param number of parameters in the current PLN model
      nb_param      = function() {(self$k-1) + sum(map_int(self$components, 'nb_param'))},
      #' @field entropy_clustering Entropy of the variational distribution of the cluster (multinomial)
      entropy_clustering = function() {-sum(.xlogx(private$tau))},
      #' #' @field entropy_latent Entropy of the variational distribution of the latent vector (Gaussian)
      #' entropy_latent       = function() {sum(self$mixtureParam * map_dbl(self$components, 'entropy'))},
      #' @field entropy Full entropy of the variational distribution (latent vector + clustering)
      entropy       = function() {self$entropy_clustering},
      #' @field loglik variational lower bound of the loglikelihood
      loglik = function() {sum(self$loglik_vec)},
      #' @field loglik_vec element-wise variational lower bound of the loglikelihood
      loglik_vec = function() {
        J_ik <- sapply(private$comp, function(comp_) comp_$loglik_vec)
        J_ik[private$tau <= .Machine$double.eps] <- 0
        rowSums(private$tau * J_ik) - rowSums(.xlogx(private$tau)) + private$tau %*% log(self$mixtureParam)
        },
      #' @field BIC variational lower bound of the BIC
      BIC        = function() {self$loglik - .5 * log(self$n) * self$nb_param},
      #' @field ICL variational lower bound of the ICL
      ICL        = function() {self$BIC - self$entropy},
      #' @field R_squared approximated goodness-of-fit criterion
      R_squared     = function() {sum(self$mixtureParam * map_dbl(self$components, "R_squared"))},
      #' @field criteria a vector with loglik, BIC, ICL, R_squared and number of parameters
      criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL, R_squared = self$R_squared)},
      #' @field model_par a list with the matrices of parameters found in the model (Theta, Sigma, plus some others depending on the variant)
      model_par  = function() {list(Theta = private$mix_up('model_par$Theta'), Sigma = private$mix_up('model_par$Sigma'))},
      #' @field var_par a list with two matrices, M and S2, which are the estimated parameters in the variational approximation
      var_par    = function() {list(M  = private$mix_up('var_par$M'), S2 = private$mix_up('var_par$S2'))},
      #' @field latent a matrix: values of the latent vector (Z in the model)
      latent = function() {private$mix_up('latent')},
      #' @field fitted a matrix: fitted values of the observations (A in the model)
      fitted = function() {private$mix_up('fitted')}
    )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR INTERNAL USE -> PLNfamily
## ----------------------------------------------------------------------
## Should only be accessed BY PLNfamily but R6 friend class don't exist

# Positions in the (euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#
# @name PLNfit_latent_pos
#
# @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
# @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
# PLNmixturefit$set("public", "latent_pos",
# function(covariates, offsets) {
#   latentPos <- private$Mu + private$M + tcrossprod(covariates, private$Theta) + offsets
#   latentPos
# })

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

