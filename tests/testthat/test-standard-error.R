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
  expected.fim.diag <- c(0.05471180694583303, 0.06872884883152279, 3.76086877733480796, 0.11120927881412997,
                         127.54002606041281354, 2.30895254021581353, 0.26961528195857087, 0.26752093401826549,
                         0.14340150907846613, 2.39389468779192738, 3.95345036877447420, 1.03589673894237277,
                         4.06712784736765087, 2.89669760144981092, 9.91782992809075026, 0.18207931687995424,
                         6.14064716477938699)

  ## Consistency of the standard error matrix
  sem <- standard_error(myPLN) %>% as.numeric()
  ## Values computed on the 2018/12/11 with PLNmodels version 0.5.9601)
  # old.expected.sem <- c(0.577407423403546, 0.577284617461014, 0.0739228099688871, 0.40821807394677,
  #                   0.0129234699024801, 0.0958134855472534, 0.267248717630853, 0.267273696185322,
  #                   0.378113801869815, 0.0928473302527288, 0.072725644559697, 0.138682400064212,
  #                   0.0723054848787022, 0.0866042221012381, 0.0461136022101119, 0.333358395876535,
  #                   0.058620515251328)
  ## Values computed with PLN dev version 0.10.0-9000 [sha: 3c28c43c80fc5ddc7b0bab60a21ddcc0058a7ba6]
  expected.sem <- c(0.61074712438904988, 0.54491959490336828, 0.07366443628152311, 0.42838223050364482,
                    0.01264965582530497, 0.09401441691780363, 0.27512477960959147, 0.27619961823650380,
                    0.37724639574619351, 0.09233140398888282, 0.07184785624702307, 0.14036011877194832,
                    0.07083665544956425, 0.08393640803448287, 0.04536215013531690, 0.33478936398005377,
                    0.05764941279130160)

  ## Temporal consistency (with previous fits of the PLN model, here fitted on the 2018/12/11 with PLNmodels version 0.5.9601)
  expect_equal(fim.diag        , expected.fim.diag, tolerance = tol)
  expect_equal(sem             , expected.sem     , tolerance = tol)

})
