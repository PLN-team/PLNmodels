context("test-standard-error")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## Error messages ---------------------
test_that("Check that fisher and standard_error return objects with proper dimensions and sign",  {

  myPLN_cov <- PLN(Abundance ~ Wind + offset(log(Offset)), data = trichoptera, control = PLN_param(config_post = list(variational_var = TRUE)))
  expect_is(myPLN_cov, "PLNfit")
  p <- myPLN_cov$p
  d <- myPLN_cov$d


  sem <- standard_error(myPLN_cov)
  ## Dimensions
  expect_equal(dim(sem), c(d, p))

  ## Names
  expect_equal(rownames(sem), rownames(coef(myPLN_cov)))
  expect_equal(colnames(sem), colnames(coef(myPLN_cov)))

  ## Standard errors are all positive
  for (i in 1:(p*d)) {
    expect_gte(sem[i], 0)
  }

})

## Fit model without covariates
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, control = PLN_param(config_post = list(variational_var = TRUE)))

## Consistency -----------------------
test_that("Check internal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-8

  n <- nrow(myPLN$fitted)
  ## Consistency of the standard error matrix
  sem <- (sqrt(n) * standard_error(myPLN)) %>% as.numeric()
  manual.sem <- 1/colMeans(myPLN$fitted) %>% sqrt()

  ## Internal consistency
  expect_equal(sem, manual.sem, tolerance = tol)

})


test_that("Check temporal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-2

  n <- nrow(myPLN$fitted)
  ## Consistency of the diagonal of the fisher matrix
  fim.diag <- 1/(n * standard_error(myPLN)^2)
  ## Values computed on the 2018/12/11 with PLNmodels version 0.5.9601)
  expected.fim.diag <- c(0.0612123698810698, 0.0612384161054906, 3.73462487824109, 0.122467107738817,
                         122.19280897578, 2.2230572191967, 0.285741065637069, 0.285687659219944,
                         0.142744327711051, 2.36736421753514, 3.85859113231971, 1.06111199011525,
                         3.90356517005791, 2.72098275756987, 9.59722821630398, 0.183645852556891,
                         5.93888146445577)

  ## Consistency of the standard error matrix
  sem <- standard_error(myPLN) %>% as.numeric()
  ## Values computed on the 2018/12/11 with PLNmodels version 0.5.9601)
  expected.sem <-
    c(0.577407423403546, 0.577284617461014, 0.0739228099688871, 0.40821807394677,
      0.0129234699024801, 0.0958134855472534, 0.267248717630853, 0.267273696185322,
      0.378113801869815, 0.0928473302527288, 0.072725644559697, 0.138682400064212,
      0.0723054848787022, 0.0866042221012381, 0.0461136022101119, 0.333358395876535,
      0.058620515251328)

  ## Temporal consistency (with previous fits of the PLN model, here fitted on the 2018/12/11 with PLNmodels version 0.5.9601)
  expect_equal(sem             , expected.sem     , tolerance = tol)

})

params <- PLNmodels:::create_parameters(n = 30, p = 3, d = 1, depths = 1e3)
B <- params$B
X <- as.matrix(params$X)
Y <- rPLN(n = nrow(X), mu = X %*% B, Sigma = params$Sigma, depths = params$depths)
log_O <- attr(Y, "offsets")
rownames(Y) <- rownames(X) <- rownames(log_O) <-1:nrow(X)
data <- prepare_data(Y, X, offset = exp(log_O))

test_that("Check that variance estimation are coherent in PLNfit",  {

  ## using postTreatments function to compute variances a posteriori
  myPLN <- PLN(Y ~ X + 0 + offset(log_O))

  config_post <-
    list(
      jackknife       = TRUE,
      bootstrap       = 50L,
      variational_var = TRUE,
      rsquared        = FALSE,
      trace           = 2
    )

  myPLN$postTreatment(Y, X, exp(log_O), config = config_post)

  tr_variational <- sum(standard_error(myPLN, "variational")^2)
  tr_bootstrap   <- sum(standard_error(myPLN, "bootstrap")^2)
  tr_jackknife   <- sum(standard_error(myPLN, "jackknife")^2)

  expect_gt(tr_variational, 0)
  expect_gt(tr_jackknife  , 0)
  expect_gt(tr_bootstrap  , 0)

  ## using control parameters
  myPLN_prime <- PLN(Abundance ~ Var_1 + 0 + offset(log(Offset)), data = data, control = PLN_param(config_post = config_post))

  tr_variational <- sum(standard_error(myPLN_prime, "variational")^2)
  tr_bootstrap   <- sum(standard_error(myPLN_prime, "bootstrap")^2)
  tr_jackknife   <- sum(standard_error(myPLN_prime, "jackknife")^2)

  expect_gt(tr_variational, 0)
  expect_gt(tr_jackknife  , 0)
  expect_gt(tr_bootstrap  , 0)
})

test_that("Check that variance estimation are coherent in PLNnetwork",  {
  myPCAs <- PLNPCA(Abundance ~ Var_1 + 0 + offset(log(Offset)), data = data, ranks = 1:3)
  myPCA <- myPCAs$models[[2]]
  B <- coef(myPCA); B[ , ] <- NA
  expect_equal(suppressWarnings(standard_error(myPCA)), B)
  expect_warning(standard_error(myPCA))
})

test_that("Check that variance estimation are coherent in PLNnetwork",  {
  myNetworks <- PLNnetwork(Abundance ~ Var_1 + 0 + offset(log(Offset)), data = data)
  myNet <- myNetworks$models[[1]]
  B <- coef(myNet); B[ , ] <- NA
  expect_equal(suppressWarnings(standard_error(myNet)), B)
  expect_warning(standard_error(myNet))
})

test_that("Check that variance estimation are coherent in PLNmixture",  {
  myModels <- PLNmixture(Abundance ~ Var_1 + 0 + offset(log(Offset)), data = data, clusters = 1:3, control = PLNmixture_param(smoothing = "none"))
  myModel <- myModels$models[[1]]
  B <- coef(myModel); B[ , ] <- NA
  expect_equal(suppressWarnings(standard_error(myModel)), B)
  expect_warning(standard_error(myModel))
})
