context("test-plnlda")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

data(mollusk)
suppressWarnings(mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate))

test_that("Check that PLN is running and robust",  {

  model1 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)

  expect_is(model1, "PLNLDAfit")
  expect_is(model1, "PLNfit")

  model2 <- PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)
  model3 <- PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)

  expect_equal(PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)$fitted,
               PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)$fitted)

  ## add more
})

test_that("Check PLNLDA initialization",  {
  tol <- 1e-4

  ## use default initialization (LM)
  model1 <- PLNLDA(Abundance ~ 0, grouping = Group, data = trichoptera)

  ## initialization with the previous fit
  model2 <- PLNLDA(Abundance ~ 0, grouping = Group,
                   data = trichoptera, control = list(inception = model1, trace = 0))

  expect_equal(model2$loglik   , model1$loglik   , tolerance = tol)
  expect_equal(model2$model_par, model1$model_par, tolerance = tol)
  expect_equal(model2$var_par  , model1$var_par  , tolerance = tol)
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
  lda_with <- PLNLDA(Abundance ~ 1 + duration + offset(log(Offset)),
                     grouping = site,
                     data = mollusc)
  lda_wo <- PLNLDA(Abundance ~ 0 + duration + offset(log(Offset)),
                     grouping = site,
                     data = mollusc)
  expect_equal(lda_with$group_means, lda_wo$group_means)
  expect_equal(coef(lda_with), coef(lda_wo))
  expect_equal(vcov(lda_with), vcov(lda_wo))
})

## add tests for predictions, tests for fit

test_that("Predictions have the right dimensions.", {
  model1 <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                   grouping = Group, data = trichoptera)

  expect_length(predict(model1, newdata = trichoptera, type = "response"),
                nrow(trichoptera))
  expect_is(predict(model1, newdata = trichoptera, type = "response"),
            "factor")
  expect_equal(dim(predict(model1, newdata = trichoptera)),
               c(nrow(trichoptera), length(levels(trichoptera$Group))))
  expect_equal(dim(predict(model1, newdata = trichoptera, type = "scores")),
               c(nrow(trichoptera), model1$rank))
  ## log-posterior probabilities are nonpositive
  expect_lt(max(predict(model1, newdata = trichoptera)), 0)
  ## Posterior probabilities are between 0 and 1
  expect_lte(max(predict(model1, newdata = trichoptera, scale = "prob")), 1)
  expect_gte(min(predict(model1, newdata = trichoptera, scale = "prob")), 0)

})


