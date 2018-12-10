context("test-plcpca")

## get data with some offset
data(trichoptera)
trichoptera$TotalCount <-
  matrix(
    rowSums(trichoptera$Abundance),
    nrow = nrow(trichoptera$Abundance),
    ncol = ncol(trichoptera$Abundance)
  )

test_that("PLNPCA runs", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera)
  expect_is(models, "PLNPCAfamily")
  ## add some...
})

