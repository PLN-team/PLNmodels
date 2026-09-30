context("test-zipln")
require(purrr)

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])

test_that("ZIPLN: Check that ZIPLN is running and robust",  {

  expect_is(zi_single <- ZIPLN(Abundance ~ 1, data = trichoptera), "ZIPLNfit")
  expect_is(zi_row <- ZIPLN(Abundance ~ 1, data = trichoptera, zi = "row"), "ZIPLNfit")
  expect_is(zi_col <- ZIPLN(Abundance ~ 1, data = trichoptera, zi = "col"), "ZIPLNfit")
  expect_is(zi_covar  <- ZIPLN(Abundance ~ 1 | Wind, data = trichoptera), "ZIPLNfit")
  expect_is(zi_covar  <- ZIPLN(Abundance ~ 0 + Wind | Wind, data = trichoptera), "ZIPLNfit")

  ## initialization could not work without any regressor...
  expect_error(ZIPLN(Abundance ~ 0, data = trichoptera))

  expect_is(ZIPLN(trichoptera$Abundance ~ 1), "ZIPLNfit")

  expect_equal(ZIPLN(trichoptera$Abundance ~ 1 + trichoptera$Wind)$fitted,
               ZIPLN(Abundance ~ Wind, data = trichoptera)$fitted, tolerance = 1e-3)

  expect_is(model_sparse <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(penalty = 0.2, trace = 0)), "ZIPLNfit_sparse")
  expect_is(model_fixed <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(Omega = zi_single$model_par$Omega, trace = 0)), "ZIPLNfit_fixed")

})

test_that("ZIPLN: Routine comparison between the different covariance models",  {
  model_full      <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "full"     , trace = 0))
  model_diagonal  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "diagonal" , trace = 0))
  model_spherical <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "spherical", trace = 0))
  expect_gte(model_full$loglik  , model_diagonal$loglik)
  expect_gte(model_diagonal$loglik, model_spherical$loglik)
})

test_that("PLN is working with a single variable data matrix",  {
  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  colnames(Y) <- "Y"
  expect_is(ZIPLN(Y ~ 1), "ZIPLNfit")

  Y <- matrix(rpois(10, exp(0.5)), ncol = 1)
  expect_is(ZIPLN(Y ~ 1), "ZIPLNfit")
})

test_that("PLN is working with unnamed data matrix",  {
  n = 15; d = 2; p = 4
  Y <- matrix(rpois(n*p, 1), n, p)
  X <- matrix(rnorm(n*d), n, d)
  expect_is(ZIPLN(Y ~ X), "ZIPLNfit")
})

 test_that("ZIPLN is working with different optimization algorithm in NLopt",  {

    MMA    <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(backend = "nlopt", config_optim = list(algorithm = "MMA")))
    CCSAQ  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(backend = "nlopt", config_optim = list(algorithm = "CCSAQ")))
    LBFGS  <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(backend = "nlopt", config_optim = list(algorithm = "LBFGS")))

    expect_equal(MMA$loglik, CCSAQ$loglik, tolerance = 1e-1) ## Almost equivalent, CCSAQ faster

    expect_error(ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(backend = "nlopt", config_optim = list(algorithm = "nawak"))))
 })


test_that("ZIPLN: Check that univariate ZIPLN models works, with matrix of numeric format",  {
  expect_no_error(uniZIPLN <- ZIPLN(Abundance[,1,drop=FALSE] ~ 1, data = trichoptera))
  expect_no_error(uniZIPLN <- ZIPLN(Abundance[,1] ~ 1, data = trichoptera))
   y <- trichoptera$Abundance[,1]
   expect_no_error(uniZIPLN <- ZIPLN(y ~ 1))
})

test_that("ZIPLN: Check that all univariate ZIPLN models are equivalent with the multivariate diagonal case",  {

  p <- ncol(trichoptera$Abundance)
  Offset <- trichoptera$Offset
  Wind <- trichoptera$Wind

  univariate_full <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(trace = 0))
  })

  univariate_diagonal <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(covariance = "diagonal", trace = 0))
  })

  univariate_spherical <- lapply(1:p, function(j) {
    Abundance <- trichoptera$Abundance[, j, drop = FALSE]
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, control = ZIPLN_param(covariance = "spherical", trace = 0))
  })

  multivariate_diagonal <-
    ZIPLN(Abundance ~ 1 + offset(log(Offset)) | Wind, data = trichoptera, control = ZIPLN_param(covariance = "diagonal", trace = 0))

  expect_true(all.equal(
    map_dbl(univariate_spherical, "nb_param"),
    map_dbl(univariate_full     , "nb_param")
  ))

  expect_true(all.equal(
    map_dbl(univariate_spherical, "nb_param"),
    map_dbl(univariate_diagonal , "nb_param")
  ))
  expect_true(all.equal(
    map_dbl(univariate_full     , "nb_param"),
    map_dbl(univariate_diagonal , "nb_param")
  ))

  expect_true(all.equal(
    map_dbl(univariate_full, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

  expect_true(all.equal(
    map_dbl(univariate_diagonal, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

   expect_true(all.equal(
    map_dbl(univariate_spherical, "loglik") %>% sum(),
    multivariate_diagonal$loglik, tolerance = 1e-2)
  )

  expect_true(all.equal(
    map(univariate_spherical, sigma) %>% map_dbl(as.double),
    map(univariate_diagonal , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

  expect_true(all.equal(
    map(univariate_spherical, sigma) %>% map_dbl(as.double),
    map(univariate_full , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

  expect_true(all.equal(
    map(univariate_diagonal, sigma) %>% map_dbl(as.double),
    map(univariate_full , sigma) %>% map_dbl(as.double), tolerance = .25
  ))

})

test_that("ZIPLN: the formula can be passed as a variable (#187)",  {
  ## full data set, so that no level of Group is empty
  utils::data("trichoptera", package = "PLNmodels", envir = environment())
  tri  <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
  ctrl <- ZIPLN_param(trace = 0)
  for (f in list(Abundance ~ 1, Abundance ~ 1 + Wind, Abundance ~ 1 + Group,
                 Abundance ~ 1 + Wind | 1 + Group)) {
    expect_is(fit <- ZIPLN(f, tri, control = ctrl), "ZIPLNfit")
    expect_equal(fit$loglik, do.call(ZIPLN, list(f, tri, control = ctrl))$loglik)
    expect_equal(dim(predict(fit, newdata = tri[1:3, ])), c(3L, ncol(tri$Abundance)))
  }
})

## Many zeros that are not excess zeros: the weaker species of each pair of
## competitors is set to 0 where the stronger one dominates (#185, #186)
simulate_competing_community <- function(seed, n = 50, p = 20, width = 10) {
  set.seed(seed)
  x     <- seq(0, 100, length.out = n)
  theta <- runif(p, 0, 100)
  mu    <- log(100) - outer(x, theta, "-")^2 / (2 * width^2)
  Y     <- matrix(rpois(n * p, exp(mu + matrix(rnorm(n * p, 0, 0.5), n, p))), n, p)
  pairs <- matrix(sample(p), ncol = 2)
  for (k in seq_len(nrow(pairs))) {
    a <- pairs[k, 1]; b <- pairs[k, 2]
    Y[Y[, a] > Y[, b], b] <- 0
  }
  keep <- colSums(Y) > 0; rows <- rowSums(Y[, keep]) > 0
  colnames(Y) <- paste0("sp", 1:p)
  suppressWarnings(prepare_data(Y[rows, keep], data.frame(x = x[rows])))
}

test_that("ZIPLN: the objective decreases and the fit is at least as good as PLN's (#185, #186)",  {
  dat <- simulate_competing_community(1)
  pln <- PLN(Abundance ~ 1 + x + I(x^2), dat, control = PLN_param(trace = 0))
  for (backend in c("builtin", "nlopt")) {
    zi <- ZIPLN(Abundance ~ 1 + x + I(x^2), dat, zi = "col",
                control = ZIPLN_param(backend = backend, trace = 0))
    obj <- zi$optim_par$objective
    expect_true(all(diff(obj) <= 1e-6 * abs(obj[-1])))
    expect_equal(zi$optim_par$objective_increases, 0L)
    expect_equal(-zi$loglik, min(obj), tolerance = 1e-6)
  }
  ## ZIPLN nests PLN (pi -> 0)
  zi <- ZIPLN(Abundance ~ 1 + x + I(x^2), dat, zi = "col", control = ZIPLN_param(trace = 0))
  expect_gte(zi$loglik, pln$loglik)
})

test_that("ZIPLN: the objective of a sparse fit includes the penalty and decreases",  {
  dat <- simulate_competing_community(1)
  zi <- ZIPLN(Abundance ~ 1 + x + I(x^2), dat, control = ZIPLN_param(penalty = 0.2, trace = 0))
  obj <- zi$optim_par$objective
  expect_true(all(diff(obj) <= 1e-6 * abs(obj[-1])))
  expect_equal(tail(obj, 1),
               -zi$loglik + .5 * zi$n * zi$penalty * sum(abs(zi$penalty_weights * zi$model_par$Omega)))
})
