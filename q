[1mdiff --git a/R/PLNblockbis.R b/R/PLNblockbis.R[m
[1mindex baec7d88..e620c04e 100644[m
[1m--- a/R/PLNblockbis.R[m
[1m+++ b/R/PLNblockbis.R[m
[36m@@ -19,7 +19,7 @@[m
 #' @examples[m
 #' data(trichoptera)[m
 #' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)[m
[31m-#' myPLN <- PLNblockbis(Abundance ~ 1, nb_blocks = 2:10, data = trichoptera, control = PLNblockbis_param(backend="nlopt-vem"))[m
[32m+[m[32m#' myPLN <- PLNblockbis(Abundance ~ 1, nb_blocks = 1:10, data = trichoptera)[m
 #' @seealso The classes [`PLNblockbisfamily`] and [`PLNblockbisfit`], and the and the configuration function [PLNblockbis_param()].[m
 #' @importFrom stats model.frame model.matrix model.response model.offset model.weights terms[m
 #' @export[m
[36m@@ -70,8 +70,9 @@[m [mPLNblockbis <- function(formula, nb_blocks = 1:5, sparsity = 0, data, subset, we[m
 PLNblockbis_param <- function([m
     backend       = c("nlopt-vem", "nlopt", "torch"),[m
     trace         = 1,[m
[31m-    init_cl       = "clustofvar",[m
[32m+[m[32m    init_cl       = "kmeans",[m
     fixed_cl      = FALSE,[m
[32m+[m[32m    route         = c("flat", "sequential"),[m
     config_optim  = list(),[m
     config_post   = list(),[m
     inception     = "lm"     # pretrained PLNfit used as initialization[m
[36m@@ -99,9 +100,9 @@[m [mPLNblockbis_param <- function([m
   }[m
 [m
   config_opt$trace <- trace[m
[31m-  config_opt$ftol_out  <- 1e-4[m
[32m+[m[32m  config_opt$ftol_out  <- 1e-5[m
   config_opt$maxit_out <- 100[m
[31m-  config_opt$route     <- "flat"[m
[32m+[m[32m  config_opt$route     <- match.arg(route)[m
   config_opt$fixed_cl  <- fixed_cl[m
   config_opt[names(config_optim)] <- config_optim[m
 [m
[1mdiff --git a/R/PLNblockbisfamily-class.R b/R/PLNblockbisfamily-class.R[m
[1mindex 9b1d696b..3d2d8cee 100644[m
[1m--- a/R/PLNblockbisfamily-class.R[m
[1m+++ b/R/PLNblockbisfamily-class.R[m
[36m@@ -57,15 +57,14 @@[m [mPLNblockbisfamily <- R6Class([m
         control_init <- control[m
         control_init$config_optim <- config_default_nlopt[m
         control_init$backend <- "nlopt"[m
[31m-        control_here <- control_init[m
[31m-        control_here$covariance <- "diagonal"[m
[31m-        myPLN_init <- PLNfit$new(responses, covariates, offsets, weights, formula, control_here)[m
[31m-        myPLN_init$optimize(responses, covariates, offsets, weights, control_here$config_optim)[m
[31m-        control$inception <- myPLN_init[m
 [m
[31m-        myPLN_init_notdiag <- PLNfit$new(responses, covariates, offsets, weights, formula, control=PLN_param(backend="nlopt", inception="lm"))[m
[31m-        myPLN_init_notdiag $optimize(responses, covariates, offsets, weights, control$config_optim)[m
[31m-        control$inceptionnotdiag <- myPLN_init_notdiag[m
[32m+[m[32m        myPLN_init_diag <- PLNfit_diagonal$new(responses, covariates, offsets, weights, formula, control_init)[m
[32m+[m[32m        myPLN_init_diag$optimize(responses, covariates, offsets, weights, control_init$config_optim)[m
[32m+[m[32m        control$inception <- myPLN_init_diag[m
[32m+[m
[32m+[m[32m        myPLN_init_full <- PLNfit$new(responses, covariates, offsets, weights, formula, control_init)[m
[32m+[m[32m        myPLN_init_full$optimize(responses, covariates, offsets, weights, control_init$config_optim)[m
[32m+[m[32m        control$inception_full <- myPLN_init_full[m
       }[m
 [m
       ## ==================================================[m
[36m@@ -73,6 +72,7 @@[m [mPLNblockbisfamily <- R6Class([m
       ##[m
       ## Either user defined or obtained (the default)[m
       ## by kmeans on the variational parameters of the means of a fully parametrized PLN[m
[32m+[m[32m      Means <- myPLN_init_full$latent_pos[m
       if (is.list(control$init_cl)) {[m
         stopifnot(length(control$init_cl) == length(nb_blocks),[m
                   all(sapply(control$init_cl, length) == private$p))[m
[36m@@ -80,19 +80,24 @@[m [mPLNblockbisfamily <- R6Class([m
       } else {[m
         blocks <- switch(control$init_cl,[m
           "kmeans" = {[m
[31m-            Means <- t(myPLN_init$var_par$M)[m
             blocks <- lapply(nb_blocks, function(k) {[m
               if (k == private$p) res <- 1:private$p[m
[31m-              else res <- kmeans(Means, centers = k, nstart = 30)$cl[m
[32m+[m[32m              else res <- kmeans(t(Means), centers = k, nstart = 30)$cl[m
[32m+[m[32m              res[m
[32m+[m[32m            })[m
[32m+[m[32m          },[m
[32m+[m[32m          "kmeansvar" = {[m
[32m+[m[32m            blocks <- lapply(nb_blocks, function(k) {[m
[32m+[m[32m              if (k == private$p) res <- 1:private$p[m
[32m+[m[32m              else res <- kmeansvar(Means, init = k, nstart = 30)$cl[m
               res[m
             })[m
           },[m
[31m-          "clustofvar" = hclustvar(myPLN_init$var_par$M) %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list(),[m
[31m-          "hclust" = hclust(as.dist(1 - cov2cor(myPLN_init$model_par$Sigma)), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list(),[m
[32m+[m[32m          "clustofvar" = hclustvar(Means) %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list(),[m
[32m+[m[32m          "hclust" = hclust(as.dist(1 - cov2cor(crossprod(Means))), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list(),[m
           "sbm" = {[m
[31m-            Sigma <- myPLN_init$model_par$Sigma[m
[32m+[m[32m            mySBM = estimateSimpleSBM(cov(Means), "gaussian", estimOption=list(verbosity=0, exploreMin=max(nb_blocks)))[m
             blocks <- lapply(nb_blocks, function(k) {[m
[31m-              mySBM = Sigma %>%  estimateSimpleSBM("gaussian", estimOption=list(verbosity=0, exploreMin=k))[m
               mySBM$setModel(k)[m
               res <-  mySBM$memberships[m
               res[m
[36m@@ -136,14 +141,12 @@[m [mPLNblockbisfamily <- R6Class([m
             kmeansvar(init = self$models[[models_order[m+1]]]$nb_block, nstart = 30) %>%[m
             pluck("cluster") %>% as_indicator() %>% .check_boundaries()[m
           self$models[[models_order[m+1]]]$update([m
[31m-            B = self$models[[models_order[m]]]$model_par$B,[m
[31m-            M = self$models[[models_order[m]]]$var_par$M %*% blocks ,[m
[31m-            S = self$models[[models_order[m]]]$var_par$S %*% blocks,[m
[31m-[m
[31m-            ###[m
[31m-            Mu = self$models[[models_order[m]]]$var_par$M,[m
[31m-            Delta = self$models[[models_order[m]]]$var_par$S,[m
[31m-[m
[32m+[m[32m            B     = self$models[[models_order[m]]]$model_par$B,[m
[32m+[m[32m            D     = self$models[[models_order[m]]]$model_par$D,[m
[32m+[m[32m            Mu    = self$models[[models_order[m]]]$var_par$Mu,[m
[32m+[m[32m            Delta = self$models[[models_order[m]]]$var_par$Delta,[m
[32m+[m[32m            M     = self$models[[models_order[m]]]$var_par$M %*% blocks,[m
[32m+[m[32m            S     = self$models[[models_order[m]]]$var_par$S %*% blocks[m
           )[m
         }[m
         if (config$trace > 1) {[m
[1mdiff --git a/R/PLNblockbisfit-class.R b/R/PLNblockbisfit-class.R[m
[1mindex 8a7aede1..e2fc24a1 100644[m
[1m--- a/R/PLNblockbisfit-class.R[m
[1m+++ b/R/PLNblockbisfit-class.R[m
[36m@@ -42,16 +42,13 @@[m [mPLNblockbisfit <- R6Class([m
 [m
       ## Initial memberships/blocks[m
       ## Overwrite PLNfit Variational parameters (dimension q)[m
[31m-      private$Delta   <- private$S[m
[31m-      private$Mu   <- private$M + covariates %*% private$B[m
[31m-      n = nrow(responses)[m
[31m-      Q = ncol(blocks)[m
[31m-[m
[31m-      private$M <- control$inceptionnotdiag$var_par$M %*% blocks[m
[31m-      private$S <- control$inceptionnotdiag$var_par$S %*% blocks[m
[31m-      # private$M   <- private$M %*% blocks[m
[31m-      # private$S   <- private$S %*% blocks[m
[32m+[m[32m      private$Delta <- control$inception$var_par$S[m
[32m+[m[32m      private$Mu    <- control$inception$var_par$M + covariates %*% control$inception$model_par$B[m
[32m+[m[32m      n <- nrow(responses)[m
[32m+[m[32m      Q <- ncol(blocks)[m
 [m
[32m+[m[32m      private$M <- control$inception$var_par$M %*% blocks[m
[32m+[m[32m      private$S <- control$inception$var_par$S %*% blocks[m
 [m
       private$Tau <- t(blocks)[m
 [m
[36m@@ -239,12 +236,20 @@[m [mPLNblockbisfit <- R6Class([m
   ),[m
   active = list([m
     ####################################[m
[31m-    #' @field blockbis_model_par to add D as a model par[m
[31m-    blockbis_model_par  = function() {list(D = private$D)},[m
[31m-    #' @field blockbis_var_par to add Mu, Delta as var par[m
[31m-    blockbis_var_par  = function() {list(Mu = private$Mu, Delta = private$Delta)},[m
[31m-    ####################################[m
 [m
[32m+[m[32m    #' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: B (covariates), Sigma (covariance), Omega (precision) and C (loadings)[m
[32m+[m[32m    model_par = function() {[m
[32m+[m[32m      par <- super$model_par[m
[32m+[m[32m      par$D <- private$D[m
[32m+[m[32m      par[m
[32m+[m[32m    },[m
[32m+[m[32m    #' @field var_par a list with the matrices of the variational parameters: M and Mu (means) and S and Delta (sqrt(variances))[m
[32m+[m[32m    var_par = function() {[m
[32m+[m[32m      par <- super$var_par[m
[32m+[m[32m      par$Mu <- private$Mu[m
[32m+[m[32m      par$Delta <- private$Delta[m
[32m+[m[32m      par[m
[32m+[m[32m    },[m
     #' @field nb_param number of parameters in the current PLN model[m
     nb_param   = function() {as.integer(self$p * self$d + .5 * self$q * (self$q + 1) + self$q - 1) + self$p},[m
     #' @field nb_block number blocks of variables (dimension of the residual covariance)[m
[36m@@ -310,7 +315,7 @@[m [mPLNblockbisfit <- R6Class([m
 #' \dontrun{[m
 #' data(trichoptera)[m
 #' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)[m
[31m-#' myPLN <- PLNblockbis(Abundance ~ 1, data = trichoptera, nb_blocks = 1:5, sparsity = 0.1)[m
[32m+[m[32m#' myPLN <- PLNblockbis(Abundance ~ 1, data = trichoptera, nb_blocks = 1:5)[m
 #' class(myPLN)[m
 #' print(myPLN)[m
 #' }[m
[1mdiff --git a/R/optimize_plnblockbis.R b/R/optimize_plnblockbis.R[m
[1mindex 8af1c9c8..0735fc6b 100644[m
[1m--- a/R/optimize_plnblockbis.R[m
[1m+++ b/R/optimize_plnblockbis.R[m
[36m@@ -34,11 +34,11 @@[m [moptimize_plnblockbis <- function(data, params, config) {[m
   objective <- Inf[m
   repeat {[m
 [m
[32m+[m[32m    # M Step[m
[32m+[m[32m    optim_B <- optim_plnblockbis_B(data, new_parameters, config)[m
[32m+[m[32m    new_parameters$B <- optim_B$B[m
     optim_Omega <- optim_plnblockbis_Omega(M = new_parameters$M, S = new_parameters$S, w = data$w)[m
     new_parameters$Omega <- optim_Omega$Omega[m
[31m-    optim_D <- optim_plnblockbis_D(X = data$X, B = new_parameters$B, Mu=new_parameters$Mu,[m
[31m-                                   Delta=new_parameters$Delta, w = data$w)[m
[31m-    new_parameters$D <- optim_D$D[m
 [m
     # VE Step[m
     optim_VE <- optim_plnblockbis_VE(data, new_parameters, config)[m
[36m@@ -46,22 +46,8 @@[m [moptimize_plnblockbis <- function(data, params, config) {[m
     new_parameters$S <- optim_VE$S[m
     new_parameters$Mu <- optim_VE$Mu[m
     new_parameters$Delta <- optim_VE$Delta[m
[31m-[m
[31m-    # Alternative[m
[31m-    # optim_VE_blocks <- optim_plnblockbis_VE_blocks(data, new_parameters, config)[m
[31m-    # new_parameters$M <- optim_VE_blocks$M[m
[31m-    # new_parameters$S <- optim_VE_blocks$S[m
[31m-    # optim_VE_species <- optim_plnblockbis_VE_species(data, new_parameters, config)[m
[31m-    # new_parameters$Mu <- optim_VE_species$Mu[m
[31m-    # new_parameters$Delta <- optim_VE_species$Delta[m
[31m-[m
[31m-    if (!config$fixed_cl) new_parameters$T <- optim_plnblockbis_Tau(data, new_parameters)[m
[31m-    #print("reached Tau optim")[m
[31m-[m
[31m-    # M Step[m
[31m-    optim_B <- optim_plnblockbis_B(data, new_parameters, config)[m
[31m-    #print("reached B optim")[m
[31m-    new_parameters$B <- optim_B$B[m
[32m+[m[32m    new_parameters$D <- optim_VE$D[m
[32m+[m[32m    new_parameters$T <- optim_VE$Tau[m
 [m
     # Going next step and assessing convergence[m
     nb_iter <- nb_iter + 1[m
[36m@@ -87,8 +73,8 @@[m [moptimize_plnblockbis <- function(data, params, config) {[m
     objective  <- new_objective[m
   }[m
   out   <- new_parameters[m
[31m-  out$A <- exp(data$O + data$X %*% out$B) * (exp(out$M + .5 * out$S**2) %*% out$T)[m
[31m-  out$Z <- data$O + data$X %*% out$B + out$M %*% out$T[m
[32m+[m[32m  out$A <- optim_VE$A[m
[32m+[m[32m  out$Z <- data$O + out$Mu + out$M %*% out$T[m
   out$Sigma <- optim_Omega$Sigma[m
   ## J'ai remplacé vloglik par loglik ici[m
   out$Ji <- plnblockbis_vloglik(data, new_parameters)[m
[1mdiff --git a/inst/case_studies/oaks_tree.R b/inst/case_studies/oaks_tree.R[m
[1mindex e2b9b76e..f9b50ded 100644[m
[1m--- a/inst/case_studies/oaks_tree.R[m
[1m+++ b/inst/case_studies/oaks_tree.R[m
[36m@@ -15,10 +15,10 @@[m [msystem.time(myPLN_diagonal <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), da[m
 system.time(myPLN_spherical <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, control = PLN_param(covariance = "spherical")))[m
 [m
 ## Blockwise covariance[m
[31m-system.time(myPLN_blocks <- PLNblock([m
[31m-            Abundance ~ 0 + tree + offset(log(Offset)), nb_blocks = seq(24, 100, by = 4), data = oaks,[m
[31m-            # sparsity = 0.25 * sqrt(seq(20, 100, by = 4)),[m
[31m-            control = PLNblock_param(inception = myPLN))[m
[32m+[m[32msystem.time(myPLN_blocks <- PLNblockbis([m
[32m+[m[32m            Abundance ~ 0 + tree + offset(log(Offset)), nb_blocks = 2:10, data = oaks,[m
[32m+[m[32m            control = PLNblockbis_param(config_optim = list(xtol_rel = 1e-10, ftol_out = 1e-6)),[m
[32m+[m[32m            sparsity = 0.25 * sqrt(2:10))[m
 )[m
 myPLN_block <- getBestModel(myPLN_blocks)[m
 [m
[1mdiff --git a/src/optim_plnblockbis.cpp b/src/optim_plnblockbis.cpp[m
[1mindex daf3e6ec..529c0d49 100644[m
[1m--- a/src/optim_plnblockbis.cpp[m
[1m+++ b/src/optim_plnblockbis.cpp[m
[36m@@ -132,11 +132,6 @@[m [mRcpp::List  optim_plnblockbis_D([m
     const arma::mat & Delta,[m
     const arma::vec & w  // (n)[m
 ) {[m
[31m-  //std::cout << "optim_plnblockbis_D" << std::endl;[m
[31m-  // const arma::mat Mu2t = Mu.t() % Mu.t() ;[m
[31m-  // const arma::mat Delta2t = Delta.t() % Delta.t() ;[m
[31m-  // arma::vec d = (Mu2t + Delta2t + XB.t() % XB.t() - 2 * (Mu.t() % XB.t())) * w;[m
[31m-[m
   const arma::mat MumXB = Mu - X * B ;[m
   arma::rowvec d = w.t() * (MumXB % MumXB + Delta % Delta) / accu(w);[m
 [m
[36m@@ -219,15 +214,15 @@[m [mRcpp::List optim_plnblockbis_VE_blocks([m
     arma::mat A_T = A2 % (A1 * T.t()) ;[m
     arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) ;[m
     arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;[m
[31m-    arma::rowvec omega = 1./diagvec(Omega).t() ;[m
[32m+[m[32m    arma::rowvec omega = diagvec(Omega).t() ;[m
 [m
     double objective =[m
       accu(w.t() * (A - Y % (O + Mu + M * T)))  // ELBO term for Poisson[m
       - 0.5 * accu(w.t() * log(S2))             // ELBO term for blocks[m
       + .5 * trace(Omega * nSigma) ;            // ...[m
 [m
[31m-      metadata.map<S_ID>(grad)     = diagmat(w) *  (S.each_row() / omega + S % A_T - pow(S, -1)) ; // question omega ici ?[m
[31m-      metadata.map<M_ID>(grad)  = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;[m
[32m+[m[32m      metadata.map<S_ID>(grad) = diagmat(w) *  (S.each_row() % omega + S % A_T - pow(S, -1)) ; // question omega ici ?[m
[32m+[m[32m      metadata.map<M_ID>(grad) = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;[m
 [m
       return objective;[m
   };[m
[36m@@ -245,9 +240,6 @@[m [mRcpp::List optim_plnblockbis_VE_blocks([m
 }[m
 [m
 [m
[31m-[m
[31m-[m
[31m-[m
 // [[Rcpp::export]][m
 Rcpp::List optim_plnblockbis_VE_species([m
     const Rcpp::List & data  , // List(Y, X, O, w)[m
[36m@@ -269,8 +261,6 @@[m [mRcpp::List optim_plnblockbis_VE_species([m
   const auto metadata = tuple_metadata(init_Mu, init_Delta);[m
   enum { Mu_ID, Delta_ID }; // Names for metadata indexes[m
 [m
[31m-[m
[31m-[m
   auto parameters = std::vector<double>(metadata.packed_size);[m
   metadata.map<Mu_ID>(parameters.data()) = init_Mu;[m
   metadata.map<Delta_ID>(parameters.data()) = init_Delta;[m
[36m@@ -278,13 +268,13 @@[m [mRcpp::List optim_plnblockbis_VE_species([m
   set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));[m
 [m
   const arma::mat XB = X * B ;[m
[32m+[m[32m  double w_bar = accu(w);[m
 [m
   // Optimize[m
[31m-  auto objective_and_grad = [&metadata, &M, &S, &Y, &X, &O, &T, &Omega, &w, &XB](const double * params, double * grad) -> double {[m
[32m+[m[32m  auto objective_and_grad = [&metadata, &M, &S, &Y, &X, &O, &T, &Omega, &w, &w_bar, &XB](const double * params, double * grad) -> double {[m
     const arma::mat Mu = metadata.map<Mu_ID>(params);[m
     const arma::mat Delta = metadata.map<Delta_ID>(params);[m
 [m
[31m-    double w_bar = accu(w);[m
     arma::mat Delta2 = Delta % Delta ;[m
     arma::mat MumXB = Mu - XB ;[m
     arma::mat S2 = S % S ;[m
[36m@@ -292,7 +282,6 @@[m [mRcpp::List optim_plnblockbis_VE_species([m
     arma::mat A2 = trunc_exp(M + .5 * S2) ;[m
     arma::mat A  = A1 % (A2 * T) ;[m
     arma::mat A_T = A2 % (A1 * T.t()) ;[m
[31m-    arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) ;[m
     arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;[m
     arma::rowvec omega = 1./diagvec(Omega).t() ;[m
 [m
[36m@@ -308,15 +297,17 @@[m [mRcpp::List optim_plnblockbis_VE_species([m
   };[m
   OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);[m
 [m
[31m-[m
   arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());[m
   arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());[m
[32m+[m[32m  arma::rowvec d = w.t() * (pow(Mu - XB,2) + Delta % Delta) / w_bar;[m
[32m+[m
   return Rcpp::List::create([m
     Rcpp::Named("status") = static_cast<int>(result.status),[m
     Rcpp::Named("iterations") = result.nb_iterations,[m
     Rcpp::Named("objective") = result.objective,[m
     Rcpp::Named("Mu") = Mu,[m
[31m-    Rcpp::Named("Delta") = Delta);[m
[32m+[m[32m    Rcpp::Named("Delta") = Delta,[m
[32m+[m[32m    Rcpp::Named("D") = arma::diagmat(d));[m
 }[m
 [m
 [m
[36m@@ -357,25 +348,33 @@[m [mRcpp::List optim_plnblockbis_VE([m
   set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(configuration["xtol_abs"]));[m
 [m
   const arma::mat XB = X * B ;[m
[32m+[m[32m  double w_bar = accu(w);[m
[32m+[m[32m  arma::vec log_pi = arma::trunc_log(mean(T,1));[m
 [m
   // Optimize[m
[31m-  auto objective_and_grad = [&metadata, &Y, &X, &O, &T, &Omega, &w, &XB](const double * params, double * grad) -> double {[m
[32m+[m[32m  auto objective_and_grad = [&metadata, &Y, &X, &O, &T, &Omega, &w, &w_bar, &XB, &log_pi](const double * params, double * grad) -> double {[m
     const arma::mat M = metadata.map<M_ID>(params);[m
     const arma::mat S = metadata.map<S_ID>(params);[m
     const arma::mat Mu = metadata.map<Mu_ID>(params);[m
     const arma::mat Delta = metadata.map<Delta_ID>(params);[m
 [m
[31m-    double w_bar = accu(w);[m
     arma::mat Delta2 = Delta % Delta ;[m
     arma::mat MumXB = Mu - XB ;[m
     arma::mat S2 = S % S ;[m
     arma::mat A1 = trunc_exp(O + Mu + .5 * Delta2) ;[m
     arma::mat A2 = trunc_exp(M + .5 * S2) ;[m
[32m+[m
[32m+[m[32m    arma::mat Tau = M.t() * diagmat(w) * Y - A2.t()* diagmat(w) * A1  ;[m
[32m+[m[32m    Tau.each_col() += log_pi ;[m
[32m+[m[32m    Tau.each_col( [](arma::vec& x){[m
[32m+[m[32m      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;[m
[32m+[m[32m    }) ;[m
[32m+[m
     arma::mat A  = A1 % (A2 * T) ;[m
     arma::mat A_T = A2 % (A1 * T.t()) ;[m
     arma::mat nSigma = (M.t() * (M.each_col() % w) + diagmat(w.t() * S2)) ;[m
     arma::rowvec d = w.t() * (MumXB % MumXB + Delta2) / w_bar;[m
[31m-    arma::rowvec omega = 1./diagvec(Omega).t() ;[m
[32m+[m[32m    arma::rowvec omega = diagvec(Omega).t() ;[m
 [m
     double objective =[m
       accu(w.t() * (A - Y % (O + Mu + M * T)))  // ELBO term for Poisson[m
[36m@@ -384,7 +383,7 @@[m [mRcpp::List optim_plnblockbis_VE([m
       - 0.5 * accu(w.t() * log(Delta2))         // ELBO term for species[m
       + 0.5 * w_bar * accu(log(d)) ;            // ...[m
 [m
[31m-    metadata.map<S_ID>(grad)     = diagmat(w) *  (S.each_row() / omega + S % A_T - pow(S, -1)) ; // question omega ici ?[m
[32m+[m[32m    metadata.map<S_ID>(grad)     = diagmat(w) *  ((S.each_row() % omega) + S % A_T - pow(S, -1)) ; // question omega ici ?[m
     metadata.map<Delta_ID>(grad) = diagmat(w) * ((Delta.each_row() / d) + Delta % A - pow(Delta, -1)) ;[m
     metadata.map<M_ID>(grad)  = diagmat(w) *  (M * Omega + A_T - Y * T.t()) ;[m
     metadata.map<Mu_ID>(grad) = diagmat(w) * ((MumXB.each_row() / d) + A - Y) ;[m
[36m@@ -397,6 +396,16 @@[m [mRcpp::List optim_plnblockbis_VE([m
     arma::mat S = metadata.copy<S_ID>(parameters.data());[m
     arma::mat Mu = metadata.copy<Mu_ID>(parameters.data());[m
     arma::mat Delta = metadata.copy<Delta_ID>(parameters.data());[m
[32m+[m[32m    arma::rowvec d = w.t() * (pow(Mu - XB,2) + Delta % Delta) / w_bar;[m
[32m+[m[32m    arma::mat A1 = trunc_exp(O + Mu + .5 * Delta%Delta) ;[m
[32m+[m[32m    arma::mat A2 = trunc_exp(M + .5 * S%S) ;[m
[32m+[m[32m    arma::mat Tau = M.t() * diagmat(w) * Y - A2.t() * diagmat(w) * A1  ;[m
[32m+[m[32m    Tau.each_col() += log_pi ;[m
[32m+[m[32m    Tau.each_col( [](arma::vec& x){[m
[32m+[m[32m      x = trunc_exp(x - max(x)) / sum(trunc_exp(x - max(x))) ;[m
[32m+[m[32m    }) ;[m
[32m+[m[32m    arma::mat A  = A1 % (A2 * T) ;[m
[32m+[m
     return Rcpp::List::create([m
       Rcpp::Named("status") = static_cast<int>(result.status),[m
       Rcpp::Named("iterations") = result.nb_iterations,[m
[36m@@ -404,5 +413,9 @@[m [mRcpp::List optim_plnblockbis_VE([m
       Rcpp::Named("M") = M,[m
       Rcpp::Named("S") = S,[m
       Rcpp::Named("Mu") = Mu,[m
[31m-      Rcpp::Named("Delta") = Delta);[m
[32m+[m[32m      Rcpp::Named("Delta") = Delta,[m
[32m+[m[32m      Rcpp::Named("Tau") = Tau,[m
[32m+[m[32m      Rcpp::Named("D") = diagmat(d),[m
[32m+[m[32m      Rcpp::Named("A") = A[m
[32m+[m[32m    );[m
 }[m
