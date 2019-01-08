context("test-plnnetworkfamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNnetwork methods", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  expect_equal(getBestModel(models), getBestModel(models, "BIC"))

  ## add some
})
