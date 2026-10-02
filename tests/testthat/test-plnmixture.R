context("test-plnmixture")

library(purrr)
data(trichoptera)
## use a subset to save some time
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

n <- nrow(trichoptera$Abundance)
p <- ncol(trichoptera$Abundance)
k <- 3

covariances <- c("spherical", "diagonal", "full")
## one collection per covariance model, without and with a covariate
models <- map(setNames(covariances, covariances), function(covariance)
  PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
             control = PLNmixture_param(covariance = covariance, smoothing = "none")))
models_cov <- map(setNames(covariances, covariances), function(covariance)
  PLNmixture(Abundance ~ 1 + Precipitation + offset(log(Offset)), clusters = 1:3, data = trichoptera,
             control = PLNmixture_param(covariance = covariance, smoothing = "none")))
mix_wo <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters = 1:3, data = trichoptera,
                     control = PLNmixture_param(smoothing = "none"))

## parameters of a component (means and covariance), and class of its covariance matrix
nb_param_component <- c(spherical = p + 1, diagonal = 2 * p, full = p + p * (p + 1) / 2)
class_Sigma <- c(spherical = "dgCMatrix", diagonal = "dgCMatrix", full = "matrix")

## the model with k components of a collection: classes, fields, methods and
## predictions, with d covariates
expect_mixture_fit <- function(family, covariance, d) {
  model <- getModel(family, k)

  expect_is(family, "PLNmixturefamily")
  expect_is(model , "PLNmixturefit")
  expect_length(model$components, k)
  for (component in model$components) {
    expect_is(component, "PLNfit")
    expect_equal(component$vcov_model, covariance)
  }

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
  expect_true(all(map_lgl(model$model_par$Sigma, inherits, class_Sigma[[covariance]])))

  ## fields and active bindings
  expect_equal(dim(model$model_par$Theta), c(d, p))
  expect_equal(dim(model$model_par$Mu), c(p, k))
  expect_true(all(map_lgl(model$model_par$Sigma, ~all.equal(dim(.x), c(p,p)))))
  expect_equal(length(model$mixtureParam), k)
  expect_equal(sum(model$loglik_vec), model$loglik)
  expect_lt(model$BIC, model$loglik)
  expect_gt(model$R_squared, 0)
  expect_equal(model$nb_param, p * d + (k - 1) + k * nb_param_component[[covariance]])

  ## S3 methods
  expect_equal(dim(fitted(model)), c(n, p))
  expect_equal(sigma(model), model$model_par$Sigma)
  if (d == 0) expect_equal(coef(model), matrix(0, 0, p))
  expect_equal(coef(model, "main")      , model$model_par$Theta)
  expect_equal(coef(model, "means")     , model$model_par$Mu)
  expect_equal(coef(model, "covariance"), model$model_par$Sigma)
  expect_equal(coef(model, "mixture")   , model$model_par$Pi)

  expect_true(inherits(plot(model, type = "pca"   , plot = FALSE), "ggplot"))
  expect_true(inherits(plot(model, type = "matrix", plot = FALSE), "ggplot"))

  ## R6 methods
  expect_true(inherits(model$plot_clustering_pca(plot = FALSE), "ggplot"))
  expect_true(inherits(model$plot_clustering_data(plot = FALSE), "ggplot"))

  ## Predictions, train = test
  predictions_response <- predict(model, newdata = trichoptera, type = "response")
  predictions_post     <- predict(model, newdata = trichoptera, "posterior")
  predictions_score    <- predict(model, newdata = trichoptera, type = "position")
  expect_length(predictions_response, n)
  expect_is(predictions_response, "factor")
  expect_equal(dim(predictions_post),  c(n, k))
  expect_equal(dim(predictions_score), c(n, p))
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predictions_post), 1)
  expect_gte(min(predictions_post), 0)

  ## Predictions, train != test
  test <- 1:nrow(trichoptera) < (nrow(trichoptera)/2)
  expect_equal(dim(predict(model, newdata = trichoptera[test, ], type = "posterior")), c(sum(test), k))
}

test_that("PLN works for abitrary cluster sequences when smoothing is requested", {
  expect_is(PLNmixture(
    Abundance ~ 1 + offset(log(Offset)), clusters = c(2, 4),
    data = trichoptera,
    control = PLNmixture_param(smoothing = "both")
  ), "PLNmixturefamily")
})

test_that("Check that PLNmixture is running and robust",  {

  for (family in models) {
    expect_is(plot(family, reverse = TRUE), "ggplot")
    expect_is(plot(family), "ggplot")
    expect_is(plot(family, type = 'diagnostic'), "ggplot")
  }

  expect_error(PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters =-2, data = trichoptera))

  expect_error(PLNmixture(Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)), data = trichoptera))

  ## the intercept is redundant with the group means
  expect_equal(sum(map_dbl(models$spherical$models, "loglik")), sum(map_dbl(mix_wo$models, "loglik")), tolerance = 1e-1)
  expect_equal(dim(predict(getModel(mix_wo, 3), newdata = trichoptera)),
               dim(predict(getModel(models$spherical, 3), newdata = trichoptera)))
})

for (covariance in covariances) {
  test_that(paste("PLNmixture fit, no covariate:", covariance, "model of the covariance is working"), {
    expect_mixture_fit(models[[covariance]], covariance, d = 0)
  })
  test_that(paste("PLNmixture fit, with covariate:", covariance, "model of the covariance is working"), {
    expect_mixture_fit(models_cov[[covariance]], covariance, d = 1)
  })
}

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
  expect_output(getModel(models$spherical, k)$show(),
                output,
                fixed = TRUE)
})

test_that("PLNmixture fit: the ICL is below the loglikelihood", {
  model <- getModel(models_cov$full, k)
  expect_lt(model$ICL, model$loglik)
})
