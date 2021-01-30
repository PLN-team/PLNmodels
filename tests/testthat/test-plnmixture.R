context("test-plnmixture")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)


test_that("Check that PLNmixture is running and robust",  {

  models <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), data = trichoptera)

  model <- getBestModel(models1)

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

# test_that("Use or not of the intercept does not change the result.",  {
#   mix_with <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
#                          control_main = list(iterates = 0))
#   mix_wo <- PLNmixture(Abundance ~ 0 + offset(log(Offset)), data = trichoptera,
#                        control_main = list(iterates = 0))
#   # expect_equivalent(mix_with, mix_wo)
#   # expect_equal(lda_with$group_means, lda_wo$group_means)
#   # expect_equal(coef(lda_with), coef(lda_wo))
#   # expect_equal(vcov(lda_with), vcov(lda_wo))
# })


