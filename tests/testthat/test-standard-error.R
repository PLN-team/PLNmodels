context("test-standard-error")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## Error messages ---------------------
test_that("Check that fisher and standard_error return objects with proper dimensions and sign",  {

  myPLN_cov <- myPLN_cov <- PLN(Abundance ~ Wind + offset(log(Offset)), data = trichoptera)
  expect_is(myPLN_cov, "PLNfit")
  p <- myPLN_cov$p
  d <- myPLN_cov$d

  fim <- vcov(myPLN_cov)
  sem <- standard_error(myPLN_cov)
  ## Dimensions
  expect_equal(dim(fim), c(p*d, p*d))
  expect_equal(dim(sem), c(p, d))

  ## Names
  expect_equal(rownames(sem), rownames(coef(myPLN_cov)))
  expect_equal(colnames(sem), colnames(coef(myPLN_cov)))

  ## Fisher is block diagonal
  expect_equal(inherits(fim, "dgCMatrix"), TRUE)

  ## Standard errors are all positive
  for (i in 1:(p*d)) {
    expect_gte(sem[i], 0)
  }

})

## Fit model without covariates
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

test_that("Fisher is deprecated", {

  expect_warning(fisher(myPLN),
                 "Deprecated: please use `vcov()` instead",
                 fixed = TRUE)
})

## Consistency -----------------------
test_that("Check internal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-8

  ## Consistency of the diagonal of the fisher matrix
  fim.diag <- Matrix::diag(vcov(myPLN))
  manual.fim.diag <- colSums(myPLN$fitted)
  ## Consistency of the standard error matrix
  sem <- standard_error(myPLN) %>% as.numeric()
  manual.sem <- 1/colSums(myPLN$fitted) %>% sqrt()

  ## Internal consistency
  expect_equal(fim.diag        , manual.fim.diag  , tolerance = tol)
  expect_equal(sem             , manual.sem       , tolerance = tol)

})


test_that("Check temporal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-3

  ## Consistency of the diagonal of the fisher matrix
  fim.diag <- Matrix::diag(vcov(myPLN)) / nrow(trichoptera)
  ## Values computed on the 2018/12/11 with PLNmodels version 0.5.9601)
  # old.expected.fim.diag <- c(0.0612123698810698, 0.0612384161054906, 3.73462487824109, 0.122467107738817,
  #                        122.19280897578, 2.2230572191967, 0.285741065637069, 0.285687659219944,
  #                        0.142744327711051, 2.36736421753514, 3.85859113231971, 1.06111199011525,
  #                        3.90356517005791, 2.72098275756987, 9.59722821630398, 0.183645852556891,
  #                        5.93888146445577)
  ## Values computed with PLN dev version 0.10.0-9000 [sha: 3c28c43c80fc5ddc7b0bab60a21ddcc0058a7ba6]
  expected.fim.diag <- c(0.291876865055442, 0.320201400524534, 3.12211437385584, 0.311716854368058,
                         122.511109484856, 2.31645762191283, 0.459742248806646, 0.380131602967321,
                         0.383615935916527, 2.04925593587078, 3.80235718830183, 1.00211905182128,
                         3.9994512862058, 2.78394125864362, 9.42584794313545, 0.380504604219339,
                         5.84154582009901)
  ## Consistency of the standard error matrix
  sem <- standard_error(myPLN) %>% as.numeric()
  ## Values computed on the 2018/12/11 with PLNmodels version 0.5.9601)
  # old.expected.sem <- c(0.577407423403546, 0.577284617461014, 0.0739228099688871, 0.40821807394677,
  #                   0.0129234699024801, 0.0958134855472534, 0.267248717630853, 0.267273696185322,
  #                   0.378113801869815, 0.0928473302527288, 0.072725644559697, 0.138682400064212,
  #                   0.0723054848787022, 0.0866042221012381, 0.0461136022101119, 0.333358395876535,
  #                   0.058620515251328)
  ## Values computed with PLN dev version 0.10.0-9000 [sha: 3c28c43c80fc5ddc7b0bab60a21ddcc0058a7ba6]
  expected.sem <- c(0.264424762823207, 0.252458702829496, 0.0808495404330058, 0.255871442486922,
                    0.0129066705181979, 0.0938619948948653, 0.210690402144234, 0.231704768388338,
                    0.230650096137834, 0.0997938686792207, 0.073261449024726, 0.142706022144864,
                    0.071433471162986, 0.0856193499666777, 0.0465309313401176, 0.23159117270696,
                    0.0591068841913978)

  ## Temporal consistency (with previous fits of the PLN model, here fitted on the 2018/12/11 with PLNmodels version 0.5.9601)
  expect_equal(fim.diag        , expected.fim.diag, tolerance = tol)
  expect_equal(sem             , expected.sem     , tolerance = tol)

})
