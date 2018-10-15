library(PLNmodels)
library(testthat)

data("trichoptera")

test_that("Check PLN initialization",  {
  tol <- 1e-4

  ## use default initialization (LM)
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0))

  ## initialization with the previous fit
  model2 <- PLN(Abundance ~ 1, data = trichoptera, control = list(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Check PLN weights",  {
  tol <- 1e-5

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Check PLN weights with spherical covariance",  {
  tol <- 1e-5

  ## no weights
  model1 <- PLN(Abundance ~ 1, data = trichoptera, covariance = "spherical", control = list(trace = 0))

  ## equivalent weigths
  model2 <- PLN(Abundance ~ 1, data = trichoptera, covariance = "spherical", weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Test different covariance models",  {
  model_full      <- PLN(Abundance ~ 1, data = trichoptera, covariance = "full", control = list(trace = 0))
  model_spherical <- PLN(Abundance ~ 1, data = trichoptera, covariance = "spherical", control = list(trace = 0))
})

## timings (weight/no weights -> 6/7% slower)
# res <- microbenchmark::microbenchmark(
#   full = PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0)),
#   ## equivalent weigths
#   spherical = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical", control = list(trace = 0)),
#   times = 20
# )
# summary(res)

res <- microbenchmark::microbenchmark(
  uw = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical",  control = list(trace = 0)),
  ## equivalent weigths
  w  = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical", weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0)),
  times = 20
)
summary(res)

#
# weights = PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = list(trace = 0))
#
# weights = PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)), control = list(trace = 0))
