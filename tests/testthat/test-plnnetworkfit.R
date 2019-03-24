context("test-plnnetworkfit")

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("PLNnetwork fit: check classes, getters and field access", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)

  X <- model.matrix(Abundance ~ 1, data = trichoptera)
  Y <- as.matrix(trichoptera$Abundance)
  n <- nrow(Y); p <- ncol(Y)
  O <- matrix(0, n, p)
  w <- rep(1, n)


  ## PLNnetworkfit
  myPLNfit <- getBestModel(models)

  ## fields and active bindings
  expect_equal(dim(myPLNfit$latent_pos(X, O)), dim(Y))
  expect_equal(dim(myPLNfit$model_par$Theta), c(ncol(Y), ncol(X)))
  expect_equal(dim(myPLNfit$model_par$Omega), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$model_par$Sigma), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$M), c(nrow(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$S), c(nrow(Y), ncol(Y)))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik)
  expect_lt(myPLNfit$BIC, myPLNfit$loglik)
  expect_lt(myPLNfit$EBIC, myPLNfit$loglik)
  expect_lt(myPLNfit$EBIC, myPLNfit$BIC)
  expect_gt(myPLNfit$R_squared, 0)
  expect_gt(myPLNfit$density, 0)
  expect_true(myPLNfit$penalty > 0)
  expect_true(is.data.frame(myPLNfit$criteria))
  ## FIXME: what about the variance parameters ? (+p??)
  expect_equal(myPLNfit$nb_param, p + myPLNfit$n_edges)

  ## S3 methods
  expect_equal(coefficients(myPLNfit), myPLNfit$model_par$Theta)
  expect_equal(dim(fitted(myPLNfit)), dim(Y))
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "main"), myPLNfit$fisher$mat)
  expect_equal(vcov(myPLNfit, "covariance"), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "covariance"), sigma(myPLNfit))
  expect_equal(dim(standard_error(myPLNfit)), dim(coefficients(myPLNfit)))
  expect_true(igraph::is.igraph(myPLNfit$plot_network(output = "igraph", plot = FALSE)))
  expect_true(inherits(myPLNfit$plot_network(output = "corrplot", plot = FALSE), "Matrix"))

})
