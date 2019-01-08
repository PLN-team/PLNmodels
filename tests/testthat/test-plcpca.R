context("test-plcpca")

## get data with some offset
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCA runs", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
  expect_is(models, "PLNPCAfamily")
  ## add some...
})

