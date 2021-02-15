context("test-plnmixture")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

mix_wt <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
                      control_main = list(iterates = 0))
mix_wo <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), data = trichoptera,
                     control_main = list(iterates = 0))
models <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

test_that("Check that PLNmixture is running and robust",  {

  model <- getBestModel(models)

  expect_is(models, "PLNmixturefamily")
  expect_is(model, "PLNmixturefit")
  expect_is(model$components[[1]], "PLNfit")

})

test_that("Invalid group throws an error",  {
  expect_error(PLNmixture(Abundance ~ 0 + offset(log(Offset)), clusters =-2, data = trichoptera))
})

test_that("Weights are unused in PLNmixture",  {
  expect_error(PLNmixture(Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)), data = trichoptera))
})


test_that("Use or not of the intercept does not change the final objective.",  {

  expect_equal(sum(map_dbl(mix_wt$models, "loglik")), sum(map_dbl(mix_wo$models, "loglik")))

})


test_that("PLNmixture fit: check classes, getters and field access",  {

  n <- nrow(trichoptera$Abundance)
  p <- ncol(trichoptera$Abundance)

  model <- getBestModel(models)

  expect_is(model, "PLNmixturefit")
  expect_equal(model$n, n)
  expect_equal(model$p, p)
  expect_equal(model$d, 0)
  expect_equal(model$rank, NULL)

  expect_null(coef(model))
  expect_true(coef(model, "mixture"))
  expect_lt(sum((model$group_means - model$model_par$Theta)^2), .Machine$double.eps)
  expect_equal(sigma(model), model$model_par$Sigma)

  expect_true(inherits(model, "numeric"))
  expect_true(inherits(model$group_means, "matrix"))
  expect_true(inherits(sigma(model), "matrix"))
#
#   ## fields and active bindings
#   expect_equal(dim(model$model_par$B), c(p, p))
#   expect_equal(dim(model$model_par$Sigma), c(p, p))
#   expect_equal(dim(model$var_par$M), c(n, p))
#   expect_equal(dim(model$var_par$S), c(n, p))
#   expect_equal(sum(model$loglik_vec), model$loglik)
#   expect_lt(model$BIC, model$loglik)
#   expect_lt(model$ICL, model$loglik)
#   expect_gt(model$R_squared, 0)
#   expect_equal(model$nb_param, p * (2 *g - 1))
#   expect_equal(dim(model$group_means), c(p, g))
#   expect_equal(dim(model$scores), c(n, model$rank))
#   expect_true(all(model$percent_var >= 0))
#   expect_equal(dim(model$corr_map), c(p, model$rank))
#
#   ## S3 methods
#   expect_equal(dim(fitted(model)), c(n, p))
#   expect_equal(sigma(model), model$model_par$Sigma)
#   expect_equal(vcov(model, "main"), model$fisher$mat)
#   expect_equal(vcov(model, "covariance"), model$model_par$Sigma)
#   expect_equal(vcov(model, "covariance"), sigma(model))
#
#   expect_true(inherits(plot(model, map = "variable", plot = FALSE), "ggplot"))
#   expect_true(inherits(plot(model, map = "individual", plot = FALSE), "ggplot"))
#   expect_true(inherits(plot(model, map = "both", plot = FALSE), "grob"))
#
#   ## R6 methods
#   expect_true(inherits(model$plot_correlation_map(plot = FALSE), "ggplot"))
#   expect_true(inherits(model$plot_individual_map(plot = FALSE), "ggplot"))
#   expect_true(inherits(model$plot_LDA(plot = FALSE), "grob"))
})

test_that("PLNmixture fit: check print message",  {

  output <- paste(
"Linear Discriminant Analysis for Poisson Lognormal distribution
==================================================================",
capture_output(print(as.data.frame(round(model$criteria, digits = 3), row.names = ""))),
"==================================================================
* Useful fields
    $model_par, $latent, $var_par, $optim_par
    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
* Useful S3 methods
    print(), coef(), sigma(), vcov(), fitted(), predict(), standard_error()
* Additional fields for LDA
    $percent_var, $corr_map, $scores, $group_means
* Additional S3 methods for LDA
    plot.PLNLDAfit(), predict.PLNLDAfit()",
sep = "\n"
)

  expect_output(model$show(),
                output,
                fixed = TRUE)
})

