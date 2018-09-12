library(PLNmodels)
library(testthat)

data("trichoptera")

test_that("Check PLN initialization",  {
  tol <- 1e-5

  ## use default initialization (LM)
  model1 <- PLN(abundance ~ 1, data = trichopetra, control = list(trace = 0))

  ## initialization with the previous fit
  model2 <- PLN(abundance ~ 1, data = trichopetra, control = list(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})
