context("test-plnlda")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("Check that PLNLDA is running and robust",  {

  model1 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)

  expect_is(model1, "PLNLDAfit")
  expect_is(model1, "PLNfit")
})

test_that("Invalid group throws an error",  {
  expect_error(PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                      grouping = group, data = trichoptera),
               "invalid grouping")
})

test_that("Check PLNLDA weights",  {
  tol <- 1e-4

  ## no weights
  model1 <- PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)

  ## equivalent weigths
  model2 <- PLNLDA(Abundance ~ 0, weights = rep(1.0, nrow(trichoptera)),
                   grouping = Group, data = trichoptera)

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
})

test_that("Use or not of the intercept does not change the result.",  {
  lda_with <- PLNLDA(Abundance ~ 1 + Wind + offset(log(Offset)),
                     grouping = Group,
                     data = trichoptera)
  lda_wo <- PLNLDA(Abundance ~ 0 + Wind + offset(log(Offset)),
                     grouping = Group,
                     data = trichoptera)
  expect_equal(lda_with$group_means, lda_wo$group_means)
  expect_equal(coef(lda_with), coef(lda_wo))
  expect_equal(vcov(lda_with), vcov(lda_wo))
})


