context("test-plcpca")

## get data with some offset
data(trichoptera)
trichoptera$TotalCount <-
  matrix(
    rowSums(trichoptera$Abundance),
    nrow = nrow(trichoptera$Abundance),
    ncol = ncol(trichoptera$Abundance)
  )

test_that("PLCPCA works", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(TotalCount)), data = trichoptera, ranks = 1:5)

  expect_equal(getBestModel(models), getBestModel(models, "ICL"))
  ## add some...
})

