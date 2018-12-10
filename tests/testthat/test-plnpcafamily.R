context("test-plnpcafamily")

data(trichoptera)
trichoptera$TotalCount <-
  matrix(
    rowSums(trichoptera$Abundance),
    nrow = nrow(trichoptera$Abundance),
    ncol = ncol(trichoptera$Abundance)
  )


test_that("PLNPCAfamily methods", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera)

  ## this is more PLNPCAfamilyfit testing
  expect_equal(getBestModel(models), getBestModel(models, "ICL"))
  ## add some...
})

