context("test-plcpca")

## get data with some offset
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
data(mollusk)
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)

test_that("PLNPCA runs", {

  models <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
  expect_is(models, "PLNPCAfamily")

  models <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = mollusc, ranks = 1:8, control_main = list(cores = 8))
  expect_is(models, "PLNPCAfamily")
    ## add some...
})

test_that("PLNPCA is fast on low ranks", {

  n <- 100
  p <- 1000
  lambda <- exp(rnorm(n * p))
  Y <- matrix(rpois(n * p, lambda), n, p)

##  profvis(models <- PLNPCA(Y ~ 1, ranks = 1:3))
  models <- PLNPCA(Y ~ 1, ranks = 1:3)
  expect_is(models, "PLNPCAfamily")
})
