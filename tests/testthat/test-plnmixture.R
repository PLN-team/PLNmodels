context("test-plnmixture")

library(purrr)
data(trichoptera)
## use a subset to save some time
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

mix_wt <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                      control = PLNmixture_param(smoothing = "none"))
mix_wo <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                     control = PLNmixture_param(smoothing = "none"))

models_sphr <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                          control = PLNmixture_param(covariance = "spherical", smoothing = "none"))

models_diag <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                          control = PLNmixture_param(covariance = "diagonal", smoothing = "none"))

models_full <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                          control = PLNmixture_param(covariance = "full", smoothing = "none"))

n <- nrow(trichoptera$Abundance)
p <- ncol(trichoptera$Abundance)
k <- 3
d <- 0

### ======================== GENERAL
test_that("PLN works for abitrary cluster sequences when smoothing is requested", {
  expect_is(PLNmixture(
    Abundance ~ 1 + offset(log(Offset)), clusters = c(2, 4),
    data = trichoptera,
    control = PLNmixture_param(smoothing = "both")
  ), "PLNmixturefamily")
})


### ============================================================================
###
### SPHERICAL, NO COVARIATE
model  <- getModel(models_sphr, k)

test_that("Check that PLNmixture is running and robust (spherical variant)",  {

  expect_is(models_sphr, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "spherical")
  expect_equal(model$components[[2]]$vcov_model, "spherical")
  expect_equal(model$components[[3]]$vcov_model, "spherical")

  expect_is(plot(models_sphr, reverse = TRUE), "ggplot")
  expect_is(plot(models_diag, reverse = TRUE), "ggplot")
  expect_is(plot(models_full, reverse = TRUE), "ggplot")
  expect_is(plot(models_sphr), "ggplot")
  expect_is(plot(models_diag), "ggplot")
  expect_is(plot(models_full), "ggplot")
  expect_is(plot(models_sphr, type = 'diagnostic'), "ggplot")
  expect_is(plot(models_diag, type = 'diagnostic'), "ggplot")
  expect_is(plot(models_full, type = 'diagnostic'), "ggplot")

  expect_error(PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters =-2, data = trichoptera))

  expect_error(PLNmixture(Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)), data = trichoptera))

  expect_equal(sum(map_dbl(mix_wt$models, "loglik")), sum(map_dbl(mix_wo$models, "loglik")), tolerance = 1e-1)

})

test_that("PLNmixture fit: check classes, getters and field access",  {

  expect_is(model, "PLNmixturefit")
  expect_equal(model$n, n)
  expect_equal(model$p, p)
  expect_equal(model$k, k)
  expect_equal(model$d, d)

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "dgCMatrix")))

  ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + k * (1 + p))

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model), matrix(0, 0, p))
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))
})

test_that("PLNmixture fit: check print message",  {
  output <- paste(
"Poisson Lognormal mixture model with 3 components and spherical covariances.",
"* Useful fields",
"    $posteriorProb, $memberships, $mixtureParam, $group_means",
"    $model_par, $latent, $latent_pos, $optim_par",
"    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria",
"    $component[[i]] (a PLNfit with associated methods and fields)",
"* Useful S3 methods",
"    print(), coef(), sigma(), fitted(), predict()",
sep="\n"
)
  expect_output(model$show(),
                output,
                fixed = TRUE)
})

test_that("Predictions have the right dimensions.", {
  predictions_response <- predict(model, newdata = trichoptera, type = "response")
  predictions_post     <- predict(model, newdata = trichoptera, "posterior")
  predictions_score    <- predict(model, newdata = trichoptera, type = "position")
  ## Train = Test
  expect_length(predictions_response, n)
  expect_is(predictions_response, "factor")
  expect_equal(dim(predictions_post),  c(n, k))
  expect_equal(dim(predictions_score), c(n, p))
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predictions_post), 1)
  expect_gte(min(predictions_post), 0)


  ## Train != Test
  test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
  expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})

test_that("Predictions are not affected by inclusion of an intercept.", {
  expect_equal(dim(predict(getModel(mix_wo, 3), newdata = trichoptera)),
               dim(predict(getModel(mix_wt, 3), newdata = trichoptera)))
})

### ============================================================================
###
### DIAGONAL, NO COVARIATE
model <- getModel(models_diag, k)
test_that("Diagonal model of the covariance is working", {

  expect_is(models_diag, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "diagonal")
  expect_equal(model$components[[2]]$vcov_model, "diagonal")
  expect_equal(model$components[[3]]$vcov_model, "diagonal")

  expect_equal(model$n, n)
  expect_equal(model$p, p)
  expect_equal(model$k, k)
  expect_equal(model$d, d)

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "dgCMatrix")))

  ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + 2 * k * p)

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model), matrix(0, 0, p))
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

   ## Prediction
   predictions_response <- predict(model, newdata = trichoptera, type = "response")
   predictions_post     <- predict(model, newdata = trichoptera, "posterior")
   predictions_score    <- predict(model, newdata = trichoptera, type = "position")
   ## Train = Test
   expect_length(predictions_response, n)
   expect_is(predictions_response, "factor")
   expect_equal(dim(predictions_post),  c(n, k))
   expect_equal(dim(predictions_score), c(n, p))
   ## Posterior probabilities are between 0 and 1
   expect_lte(max(predictions_post), 1)
   expect_gte(min(predictions_post), 0)

   ## Train != Test
   test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
   expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})

### ============================================================================
###
### FULL, NO COVARIATE
model <- getModel(models_full, k)
test_that("Full model of the covariance is working", {

  expect_is(models_diag, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "full")
  expect_equal(model$components[[2]]$vcov_model, "full")
  expect_equal(model$components[[3]]$vcov_model, "full")

  expect_equal(model$n, n)
  expect_equal(model$p, p)
  expect_equal(model$k, k)
  expect_equal(model$d, d)

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "matrix")))

  ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + k * (p + p * (p + 1) / 2) )

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model), matrix(0, 0, p))
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

   ## Prediction
   predictions_response <- predict(model, newdata = trichoptera, type = "response")
   predictions_post     <- predict(model, newdata = trichoptera, "posterior")
   predictions_score    <- predict(model, newdata = trichoptera, type = "position")
   ## Train = Test
   expect_length(predictions_response, n)
   expect_is(predictions_response, "factor")
   expect_equal(dim(predictions_post),  c(n, k))
   expect_equal(dim(predictions_score), c(n, p))
   ## Posterior probabilities are between 0 and 1
   expect_lte(max(predictions_post), 1)
   expect_gte(min(predictions_post), 0)

   ## Train != Test
   test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
   expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})

#### With COVARIATE

models_sphr_cov <- PLNmixture(Abundance ~ 1 + Precipitation + offset(log(Offset)),
                              clusters = 1:3, data = trichoptera,
                              control = PLNmixture_param(covariance = "spherical", smoothing = "none"))

models_diag_cov <- PLNmixture(Abundance ~ 1 + Precipitation + offset(log(Offset)),
                              clusters = 1:3, data = trichoptera,
                              control = PLNmixture_param(covariance = "diagonal", smoothing = "none"))

models_full_cov <- PLNmixture(Abundance ~ 1 + Precipitation + offset(log(Offset)),
                              clusters = 1:3, data = trichoptera,
                              control = PLNmixture_param(covariance = "full", smoothing = "none"))
d <- 1

model  <- getModel(models_sphr_cov, k)
test_that("Spherical model of the covariance is working with covariate", {

  expect_is(models_sphr, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "spherical")
  expect_equal(model$components[[2]]$vcov_model, "spherical")
  expect_equal(model$components[[3]]$vcov_model, "spherical")

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "dgCMatrix")))

   ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + k * (p + 1) )

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

   ## Prediction
   predictions_response <- predict(model, newdata = trichoptera, type = "response")
   predictions_post     <- predict(model, newdata = trichoptera, "posterior")
   predictions_score    <- predict(model, newdata = trichoptera, type = "position")
   ## Train = Test
   expect_length(predictions_response, n)
   expect_is(predictions_response, "factor")
   expect_equal(dim(predictions_post),  c(n, k))
   expect_equal(dim(predictions_score), c(n, p))
   ## Posterior probabilities are between 0 and 1
   expect_lte(max(predictions_post), 1)
   expect_gte(min(predictions_post), 0)

   ## Train != Test
   test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
   expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})

model  <- getModel(models_diag_cov, k)
test_that("Diagonal model of the covariance is working with covariate", {

  expect_is(models_sphr, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "diagonal")
  expect_equal(model$components[[2]]$vcov_model, "diagonal")
  expect_equal(model$components[[3]]$vcov_model, "diagonal")

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "dgCMatrix")))

   ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + k * (p + p) )

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

   ## Prediction
   predictions_response <- predict(model, newdata = trichoptera, type = "response")
   predictions_post     <- predict(model, newdata = trichoptera, "posterior")
   predictions_score    <- predict(model, newdata = trichoptera, type = "position")
   ## Train = Test
   expect_length(predictions_response, n)
   expect_is(predictions_response, "factor")
   expect_equal(dim(predictions_post),  c(n, k))
   expect_equal(dim(predictions_score), c(n, p))
   ## Posterior probabilities are between 0 and 1
   expect_lte(max(predictions_post), 1)
   expect_gte(min(predictions_post), 0)

   ## Train != Test
   test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
   expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})

model  <- getModel(models_full_cov, k)
test_that("Full model of the covariance is working with covariate", {

  expect_is(models_sphr, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")
  expect_is(model$components[[2]], "PLNfit")
  expect_is(model$components[[3]], "PLNfit")
  expect_equal(model$components[[1]]$vcov_model, "full")
  expect_equal(model$components[[2]]$vcov_model, "full")
  expect_equal(model$components[[3]]$vcov_model, "full")

  expect_equal(model$model_par$Mu, model$group_means)
  expect_equal(model$model_par$Sigma, sigma(model))
  expect_equal(model$model_par$Pi, model$mixtureParam)
  expect_true(inherits(model$mixtureParam   , "numeric"))
  expect_true(inherits(model$group_means    , "data.frame"))
  expect_true(inherits(model$model_par$Sigma, "list"))
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, "matrix")))

   ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(d, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_lt(model$ICL, model$loglik)
   expect_gt(model$R_squared, 0)
   expect_equal(model$nb_param, p * d + (k - 1) + k * (p + p * (p+1)/2) )

   ## S3 methods
   expect_equal(dim(fitted(model)), c(n, p))
   expect_equal(sigma(model), model$model_par$Sigma)
   expect_equal(coef(model, "main")      , model$model_par$Theta)
   expect_equal(coef(model, "means")     , model$model_par$Mu)
   expect_equal(coef(model, "covariance"), model$model_par$Sigma)
   expect_equal(coef(model, "mixture")   , model$model_par$Pi)

   expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
   expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

   ## R6 methods
   expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
   expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

   ## Prediction
   predictions_response <- predict(model, newdata = trichoptera, type = "response")
   predictions_post     <- predict(model, newdata = trichoptera, "posterior")
   predictions_score    <- predict(model, newdata = trichoptera, type = "position")
   ## Train = Test
   expect_length(predictions_response, n)
   expect_is(predictions_response, "factor")
   expect_equal(dim(predictions_post),  c(n, k))
   expect_equal(dim(predictions_score), c(n, p))
   ## Posterior probabilities are between 0 and 1
   expect_lte(max(predictions_post), 1)
   expect_gte(min(predictions_post), 0)

   ## Train != Test
   test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
   expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))

})
