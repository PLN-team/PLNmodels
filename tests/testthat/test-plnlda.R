context("test-plnlda")

data("trichoptera")

test_that("Check PLNLDA initialization",  {
  tol <- 1e-5

  ## use default initialization (LM)
  model1 <- PLNLDA(trichoptera$Abundance ~ 0, grouping = trichoptera$Group)

  ## initialization with the previous fit
  model2 <- PLNLDA(trichoptera$Abundance ~ 0, grouping = trichoptera$Group, control = list(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Check PLNLDA weights",  {
  tol <- 1e-5

  ## no weights
  model1 <- PLNLDA(trichoptera$Abundance ~ 0, grouping = trichoptera$Group)

  ## equivalent weigths
  model2 <- PLNLDA(trichoptera$Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)), grouping = trichoptera$Group)

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})
