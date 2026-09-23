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
  expect_equal(dim(myPLNfit$latent), dim(Y))
  expect_equal(dim(myPLNfit$model_par$B), c(ncol(X), ncol(Y)))
  expect_equal(dim(myPLNfit$model_par$Omega), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$model_par$Sigma), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$M), c(nrow(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$S), c(nrow(Y), ncol(Y)))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik)
  expect_lt(myPLNfit$BIC, myPLNfit$loglik)
  expect_lt(myPLNfit$EBIC, myPLNfit$loglik)
  expect_lt(myPLNfit$EBIC, myPLNfit$BIC)
  expect_equal(myPLNfit$R_squared, NA)
  expect_gt(myPLNfit$density, 0)
  expect_true(myPLNfit$penalty > 0)
  expect_true(is.data.frame(myPLNfit$criteria))
  expect_equal(myPLNfit$nb_param, 2 *p + myPLNfit$n_edges)

  ## S3 methods
  expect_equal(coefficients(myPLNfit), myPLNfit$model_par$B)
  expect_equal(dim(fitted(myPLNfit)), dim(Y))
  expect_equal(sigma(myPLNfit), myPLNfit$model_par$Sigma)
  expect_error(vcov(myPLNfit, "main"))
  expect_null(myPLNfit$vcov_coef)
  expect_equal(vcov(myPLNfit, "covariance"), myPLNfit$model_par$Sigma)
  expect_equal(vcov(myPLNfit, "covariance"), sigma(myPLNfit))
  expect_warning(standard_error(myPLNfit))
  expect_true(igraph::is_igraph(myPLNfit$plot_network(output = "igraph", plot = FALSE)))
  expect_true(inherits(myPLNfit$plot_network(output = "corrplot", plot = FALSE), "Matrix"))

})

test_that("PLNnetwork fit accepts torch backend", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())

  data("trichoptera", package = "PLNmodels", envir = environment())
  trichoptera_small <- prepare_data(
    trichoptera$Abundance[1:10, 1:4],
    trichoptera$Covariate[1:10, , drop = FALSE]
  )
  Y <- as.matrix(trichoptera_small$Abundance)
  torch_control <- PLNnetwork_param(
    backend = "torch",
    trace = 0,
    config_optim = list(
      algorithm = "RPROP",
      lr = 0.01,
      num_epoch = 5,
      num_batch = 1,
      maxit_em = 2
    )
  )

  models <- NULL
  expect_no_error(models <- PLNnetwork(
    Abundance ~ 1,
    data = trichoptera_small,
    penalties = 0.1,
    control = torch_control
  ))
  expect_false(is.null(models))

  myPLNfit <- getBestModel(models)
  expect_equal(dim(myPLNfit$latent), dim(Y))
  expect_equal(dim(myPLNfit$model_par$B), c(1, ncol(Y)))
  expect_equal(dim(myPLNfit$model_par$Omega), c(ncol(Y), ncol(Y)))
  expect_equal(dim(myPLNfit$var_par$M), dim(Y))
  expect_equal(dim(myPLNfit$var_par$S), dim(Y))
  expect_equal(sum(myPLNfit$loglik_vec), myPLNfit$loglik, tolerance = 1e-4)
})

test_that("PLNnetwork fit: graphical Lasso convergence is monitored", {
  models <- PLNnetwork(Abundance ~ 1, data = trichoptera,
                       control = PLNnetwork_param(trace = 0, n_penalties = 5))
  nonconv <- vapply(models$models, function(m) m$optim_par$glasso_nonconverged, integer(1))
  expect_equal(nonconv, rep(0L, length(models$models)))
})

test_that("the EBIC is the one of Foygel and Drton, with a tunable gamma", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)
  fit    <- getBestModel(models, "BIC")
  p <- fit$p; E <- fit$n_edges
  expect_gt(E, 0)

  ## BIC minus 2 gamma |E| log(p), the additional penalty on the edge set
  expect_equal(fit$EBIC, fit$BIC - 2 * 0.5 * E * log(p))
  expect_equal(fit$ebic_gamma, 0.5)

  ## gamma = 0 gives back the BIC, gamma = 1 penalizes twice as much as the default
  fit$ebic_gamma <- 0
  expect_equal(fit$EBIC, fit$BIC)
  fit$ebic_gamma <- 1
  expect_equal(fit$EBIC, fit$BIC - 2 * E * log(p))
  fit$ebic_gamma <- 0.5

  expect_error(fit$ebic_gamma <- 2)
  expect_error(fit$ebic_gamma <- -1)
  expect_error(fit$ebic_gamma <- NA)

  ## the collection shares a single gamma, and it drives the selection
  expect_equal(models$ebic_gamma, 0.5)
  models$ebic_gamma <- 0
  expect_equal(models$criteria$EBIC, models$criteria$BIC)
  expect_equal(getBestModel(models, "EBIC")$penalty, getBestModel(models, "BIC")$penalty)
  models$ebic_gamma <- 0.5
  ## a stronger gamma never selects a denser network than a weaker one
  sparse_at <- function(g) {models$ebic_gamma <- g; getBestModel(models, "EBIC")$n_edges}
  expect_true(all(diff(sapply(c(1, 0.5, 0), sparse_at)) >= 0))
})

test_that("the density is the proportion of edges among the possible ones", {

  models <- PLNnetwork(Abundance ~ 1, data = trichoptera)

  for (fit in models$models[c(1, 15, 30)]) {
    expect_equal(fit$density, fit$n_edges / (fit$p * (fit$p - 1) / 2))
    ## the diagonal is not a possible edge, and does not count in either term
    expect_equal(fit$density, mean(fit$latent_network("support")[upper.tri(diag(fit$p))] != 0))
  }
  expect_lte(max(models$criteria$density), 1)
})

test_that("igraph edges carry the strength of the partial correlation in their opacity", {

  fit <- getBestModel(PLNnetwork(Abundance ~ 1, data = trichoptera))
  expect_gt(fit$n_edges, 2)

  G <- fit$plot_network(output = "igraph", plot = FALSE)
  alpha <- strtoi(substr(igraph::E(G)$color, 8, 9), base = 16L) / 255
  weight <- abs(igraph::E(G)$weight)

  ## opacity increases with the strength of the edge, the strongest being opaque
  ## (only up to ties: the alpha channel is quantized to 8 bits)
  expect_false(is.unsorted(alpha[order(weight)]))
  expect_equal(max(alpha), 1, tolerance = 1e-2)
  expect_gte(min(alpha), 0.2 - 1e-8)

  ## it can be switched off, and the support has no strength to display
  G <- fit$plot_network(output = "igraph", edge.alpha = 1, plot = FALSE)
  expect_equal(unique(nchar(igraph::E(G)$color)), 9L)
  expect_equal(unique(substr(igraph::E(G)$color, 8, 9)), "FF")
  G <- fit$plot_network(type = "support", output = "igraph", plot = FALSE)
  expect_equal(unique(nchar(igraph::E(G)$color)), 7L)
})
