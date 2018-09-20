library(PLNmodels)
library(testthat)

data("trichoptera")

test_that("Check PLN initialization",  {
  tol <- 1e-5

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
  model2 <- PLN(Abundance ~ 1, data = trichoptera, weights = rep(1, nrow(trichoptera)), control = list(trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

# ## no weights
# res <- microbenchmark::microbenchmark(
#   noweights = PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0)),
#
#   ## equivalent weigths
#   weights = PLN(Abundance ~ 1, data = trichoptera, weights = rep(1, nrow(trichoptera)), control = list(trace = 0))
# )
#
weights = PLN(Abundance ~ 1, data = trichoptera, weights = runif(nrow(trichoptera)) , control = list(trace = 0))
