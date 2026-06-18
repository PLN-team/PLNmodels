context("test-pln-genpop")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
p <- ncol(trichoptera$Abundance)

## fixed p x p correlation matrix (e.g. a genetic relationship matrix), AR1-like, unit diagonal
C <- 0.5^abs(outer(1:p, 1:p, "-")); diag(C) <- 1

test_that("PLN with genpop covariance: requires C to be specified",  {
  expect_error(PLN_param(covariance = "genpop"))
  expect_error(PLN_param(covariance = "genpop", C = "not_a_matrix"))
  expect_no_error(PLN_param(covariance = "genpop", C = C))
})

test_that("PLN with genpop covariance: basic fit, builtin backend",  {
  model <- PLN(Abundance ~ 1, data = trichoptera,
               control = PLN_param(covariance = "genpop", C = C, backend = "builtin", trace = 0))

  expect_is(model, "PLNfit")
  expect_equal(model$vcov_model, "genpop")
  expect_equal(model$nb_param, p * 1 + 2)
  expect_true(is.finite(model$loglik))
  expect_equal(dim(model$model_par$Sigma), c(p, p))
  expect_equal(dim(model$model_par$Omega), c(p, p))
})

test_that("PLN with genpop covariance: basic fit, nlopt backend",  {
  model <- PLN(Abundance ~ 1, data = trichoptera,
               control = PLN_param(covariance = "genpop", C = C, backend = "nlopt", trace = 0))

  expect_is(model, "PLNfit")
  expect_equal(model$vcov_model, "genpop")
  expect_true(is.finite(model$loglik))
})

test_that("PLN with genpop covariance: builtin and nlopt backends agree",  {
  model_builtin <- PLN(Abundance ~ 1, data = trichoptera,
                        control = PLN_param(covariance = "genpop", C = C, backend = "builtin", trace = 0))
  model_nlopt   <- PLN(Abundance ~ 1, data = trichoptera,
                        control = PLN_param(covariance = "genpop", C = C, backend = "nlopt", trace = 0))

  expect_equal(model_builtin$loglik, model_nlopt$loglik, tolerance = 1e-1)
  ## small dataset (n=49): (sigma2, rho) agree loosely between backends, consistent with
  ## the loglik tolerance above -- not meant to be a tight numerical equivalence check.
  expect_equal(model_builtin$gen_par$sigma2, model_nlopt$gen_par$sigma2, tolerance = 0.2)
  expect_equal(model_builtin$gen_par$rho   , model_nlopt$gen_par$rho   , tolerance = 0.2)
})

test_that("PLN with genpop covariance: gen_par exposes valid (sigma2, rho)",  {
  model <- PLN(Abundance ~ 1, data = trichoptera,
               control = PLN_param(covariance = "genpop", C = C, trace = 0))

  expect_setequal(names(model$gen_par), c("sigma2", "rho"))
  expect_gt(model$gen_par$sigma2, 0)
  expect_gte(model$gen_par$rho, 0)
  expect_lte(model$gen_par$rho, 1)
})

test_that("PLN with genpop covariance: more parsimonious than full covariance",  {
  model_genpop <- PLN(Abundance ~ 1, data = trichoptera,
                       control = PLN_param(covariance = "genpop", C = C, trace = 0))
  model_full   <- PLN(Abundance ~ 1, data = trichoptera,
                       control = PLN_param(covariance = "full", trace = 0))

  ## genpop's Sigma is constrained to a 2-parameter family, a strict subset of all
  ## SPD matrices reachable by full covariance -- its ELBO cannot exceed full's.
  expect_lte(model_genpop$loglik, model_full$loglik)
  expect_lt(model_genpop$nb_param, model_full$nb_param)
})

test_that("PLN with genpop covariance: predict works at level 0 and level 1",  {
  model <- PLN(Abundance ~ 1, data = trichoptera,
               control = PLN_param(covariance = "genpop", C = C, trace = 0))

  newdata <- trichoptera[1:5, ]
  pred0 <- predict(model, newdata = newdata, type = "response")
  expect_equal(dim(pred0), c(5, p))
  expect_false(anyNA(pred0))

  pred1 <- predict(model, newdata = newdata, responses = trichoptera$Abundance[1:5, ], type = "response")
  expect_equal(dim(pred1), c(5, p))
  expect_false(anyNA(pred1))
})
