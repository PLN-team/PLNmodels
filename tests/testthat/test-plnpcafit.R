context("test-plnpcafit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNPCA fit: check classes, getters and field access", {
  models <- PLNPCA(Abundance ~ 1, data = trichoptera)

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  n <- nrow(Y); p <- ncol(Y)
  O <- matrix(0, nrow = n, ncol = p)
  w <- rep(1, n)

  myPLNfit <- getBestModel(models)

  ## fields and active bindings
  expect_equal(dim(myPLNfit$latent_pos(X, O)), dim(Y))
  expect_equal(dim(myPLNfit$model_par$Theta), c(ncol(Y), ncol(X)))
  expect_equal(dim(myPLNfit$model_par$B), c(ncol(Y), myPLNfit$rank))
  expect_equal(dim(myPLNfit$model_par$Sigma), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$M), c(nrow(Y), myPLNfit$rank))
  expect_equal(dim(myPLNfit$var_par$S), c(nrow(Y), myPLNfit$rank))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik)
  expect_lt(myPLNfit$BIC, myPLNfit$loglik)
  expect_lt(myPLNfit$ICL, myPLNfit$loglik)
  expect_lt(myPLNfit$ICL, myPLNfit$BIC)
  expect_gt(myPLNfit$R_squared, 0)
  expect_equal(myPLNfit$nb_param, p + p * myPLNfit$rank)
  expect_equal(dim(myPLNfit$rotation), c(p, myPLNfit$rank))
  expect_equal(dim(myPLNfit$scores), c(n, myPLNfit$rank))
  expect_true(all(myPLNfit$percent_var >= 0))
  expect_equal(dim(myPLNfit$corr_circle), c(p, myPLNfit$rank))

  ## S3 methods
  expect_equal(coefficients(myPLNfit), myPLNfit$model_par$Theta)
  expect_equal(dim(fitted(myPLNfit)), dim(Y))
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "main"), myPLNfit$fisher$mat)
  expect_equal(vcov(myPLNfit, "covariance"), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "covariance"), sigma(myPLNfit))
  expect_equal(dim(standard_error(myPLNfit)), dim(coefficients(myPLNfit)))

  expect_true(inherits(plot(myPLNfit, map = "variable", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "individual", plot = FALSE), "ggplot"))
  expect_true(inherits(plot(myPLNfit, map = "both", plot = FALSE), "grob"))

  ## R6 methods
  expect_true(inherits(myPLNfit$plot_correlation_circle(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_individual_map(plot = FALSE), "ggplot"))
  expect_true(inherits(myPLNfit$plot_PCA(plot = FALSE), "grob"))

})
