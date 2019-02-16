context("test-plnpcafamily")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCAfamily methods", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

  ## this is more PLNPCAfamilyfit testing
  expect_equal(getBestModel(models), getBestModel(models, "ICL"))
  ## add some...
})

