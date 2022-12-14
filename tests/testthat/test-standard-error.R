context("test-standard-error")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## Error messages ---------------------
test_that("Check that fisher and standard_error return objects with proper dimensions and sign",  {

  myPLN_cov <- myPLN_cov <- PLN(Abundance ~ Wind + offset(log(Offset)), data = trichoptera)
  expect_is(myPLN_cov, "PLNfit")
  p <- myPLN_cov$p
  d <- myPLN_cov$d


  sem <- standard_error(myPLN_cov)
  ## Dimensions
  expect_equal(dim(sem), c(p, d))

  ## Names
  expect_equal(rownames(sem), rownames(coef(myPLN_cov)))
  expect_equal(colnames(sem), colnames(coef(myPLN_cov)))

  ## Standard errors are all positive
  for (i in 1:(p*d)) {
    expect_gte(sem[i], 0)
  }

})

## Fit model without covariates
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

## Consistency -----------------------
test_that("Check internal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-8

  ## Consistency of the standard error matrix
  sem <- standard_error(myPLN) %>% as.numeric()
  manual.sem <- 1/colSums(myPLN$fitted) %>% sqrt()

  ## Internal consistency
  expect_equal(sem, manual.sem, tolerance = tol)

})


test_that("Check temporal consistency of Fisher matrix for PLN models with no covariates",  {
  tol <- 1e-2

  ## Consistency of the diagonal of the fisher matrix
  fim.diag <- 1/(standard_error(myPLN)^2)
  ## Values computed on the 2022/12/14 with PLNmodels version 0.11.)
  expected.fim.diag <-
    c(2.9998349542779569887, 3.0010819594189506176, 183.0098806185363571331, 6.0003428448505022885, 5987.6778666489190072753,
      108.9764130656136984499, 13.9999160967500699826, 13.9996217322142406658, 6.9990426950361328551, 116.0021002230054705251,
      189.1877628339292698456, 51.9964910661757357957, 191.1007171515619802449, 133.1215538865650103162, 470.6676096783645562027,
      8.9981115680731615925, 291.0168391344596443560)
  expect_equivalent(fim.diag   , expected.fim.diag, tolerance = tol)

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
