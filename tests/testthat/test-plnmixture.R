context("test-plnmixture")

library(purrr)
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

mix_wt <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
                      control_main = list(iterates = 0))
mix_wo <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), data = trichoptera,
                     control_main = list(iterates = 0))
models <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
n <- nrow(trichoptera$Abundance)
p <- ncol(trichoptera$Abundance)
k <- 3
d <- 0

model  <- getModel(models, k)

test_that("Check that PLNmixture is running and robust",  {

  expect_is(models, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")

  expect_error(PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters =-2, data = trichoptera))

  expect_error(PLNmixture(Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)), data = trichoptera))

  expect_equal(sum(map_dbl(mix_wt$models, "loglik")), sum(map_dbl(mix_wo$models, "loglik")))

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
  expect_true(inherits(model$model_par$Sigma, "matrix"))

  ## fields and active bindings
   expect_equal(dim(model$model_par$Theta), c(0, p))
   expect_equal(dim(model$model_par$Mu), c(p, k))
   expect_equal(dim(model$model_par$Sigma), c(p, p))
   expect_equal(length(model$mixtureParam), k)
   expect_equal(dim(model$var_par$M), c(n, p))
   expect_equal(dim(model$var_par$S), c(n, 1))
   expect_equal(sum(model$loglik_vec), model$loglik)
   expect_lt(model$BIC, model$loglik)
   expect_lt(model$ICL, model$loglik)
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
"Poisson Lognormal mixture model with 3 components.",
"* check fields $posteriorProb, $memberships, $mixtureParam and $components",
"* check methods $plot_clustering_data, $plot_clustering_pca",
"* each $component[[i]] is a PLNfit with associated methods and fields",
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
               dim(predict(getModel(mix_wt, 3) , newdata = trichoptera)))
})
