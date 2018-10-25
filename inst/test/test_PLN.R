library(PLNmodels)
library(testthat)

data("trichoptera")

test_that("Check PLN initialization - fully parametrized covariance",  {
  tol <- 1e-4

  ## use default initialization (LM)
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0))

  ## initialization with the previous fit
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = list(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Check PLN initialization - diagonal parametrized covariance",  {
  tol <- 1e-4

  ## use default initialization (LM)
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0, covariance = "diagonal"))

  ## initialization with the previous fit
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = list(inception = model1, trace = 0, covariance = "diagonal"))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Check PLN weights",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("Check PLN weights with spherical covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "spherical", trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = list(covariance = "spherical", trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("Check PLN weights with diagonal covariance",  {
  tol <- 1e-2

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "diagonal", trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = list(covariance = "diagonal", trace = 0))

  model3 <- PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = list(covariance = "diagonal", trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
})

test_that("Test different covariance models",  {
  model_full      <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "full"     , trace = 0))
  model_spherical <- PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "spherical", trace = 0))
})

## timings (full, spherical, diagonal)
res <- microbenchmark::microbenchmark(
  spherical = PLN(Abundance ~ 1, data = trichoptera, control = list(covariance  ="spherical", trace = 0)),
  diagonal  = PLN(Abundance ~ 1, data = trichoptera, control = list(covariance  ="diagonal", trace = 0)),
  full      = PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0)),
  times = 20
)
summary(res)

# res <- microbenchmark::microbenchmark(
#   uw = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical",  control = list(trace = 0)),
#   ## equivalent weigths
#   w  = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical", weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0)),
#   times = 20
# )
# summary(res)

#
# weights = PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = list(trace = 0))
#
# weights = PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = list(trace = 0))
