###############################################################################
## Scale of the l1 penalty and floor on the variational means in network fits
## (penalty_scale and latent_floor in PLNnetwork_param()).
##
## The entries of a precision matrix are not scale invariant, so that a species
## with a large latent variance has nearly free edges under a penalty on the
## covariance scale: a species often absent but abundant when present ends up
## connected to most of the others. The penalty on the correlation scale
## removes this artefact; the floor keeps the latent variances from diverging.
###############################################################################

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

test_that("the defaults are the covariance scale and no floor, and are unchanged", {
  ctrl <- PLNnetwork_param(trace = 0, n_penalties = 5)
  expect_equal(ctrl$penalty_scale, "covariance")
  expect_null(ctrl$latent_floor)
  nets <- PLNnetwork(Abundance ~ 1, trichoptera, control = ctrl)
  expect_equal(nets$models[[1]]$penalty_scale, "covariance")
  expect_null(nets$models[[1]]$latent_floor)
  explicit <- PLNnetwork(Abundance ~ 1, trichoptera,
                         control = PLNnetwork_param(trace = 0, n_penalties = 5, penalty_scale = "covariance"))
  expect_identical(nets$criteria$loglik, explicit$criteria$loglik)

  expect_error(PLNnetwork_param(latent_floor = -1), "latent_floor")
  expect_error(PLNnetwork_param(latent_floor = c(1, 2)), "latent_floor")
  expect_error(PLNnetwork_param(penalty_scale = "nawak"))
})

test_that("the penalty on the correlation scale is the graphical Lasso on the correlation matrix", {
  set.seed(1)
  p <- 8
  S <- crossprod(matrix(rnorm(30 * p), 30, p)) / 30 * tcrossprod(10^runif(p, -1, 1))
  W <- matrix(1, p, p); diag(W) <- 0
  rho <- PLNmodels:::glasso_penalty(0.2, W, S, "correlation")
  expect_equal(rho, 0.2 * W * tcrossprod(sqrt(diag(S))))
  expect_equal(PLNmodels:::glasso_penalty(0.2, W, S, "covariance"), 0.2 * W)
  ## Omega = D^-1/2 Theta D^-1/2, with Theta the solution on the correlation matrix
  d <- sqrt(diag(S))
  expect_equal(graphical_lasso(S, rho)$wi, graphical_lasso(cov2cor(S), 0.2 * W)$wi / tcrossprod(d),
               tolerance = 1e-8)
})

test_that("PLNnetwork: on the correlation scale the penalties are between 0 and 1, from the empty network", {
  nets <- PLNnetwork(Abundance ~ 1, trichoptera,
                     control = PLNnetwork_param(trace = 0, penalty_scale = "correlation"))
  expect_true(all(nets$penalties > 0 & nets$penalties <= 1))
  expect_true(all(vapply(nets$models, function(m) m$penalty_scale, character(1)) == "correlation"))
  edges <- nets$criteria$n_edges[order(nets$criteria$param, decreasing = TRUE)]
  expect_lte(edges[1], 1)
  expect_gt(max(edges), 10)
  expect_true(all(is.finite(nets$criteria$loglik)))
  expect_is(getBestModel(nets, "BIC"), "PLNnetworkfit")
})

test_that("PLNnetwork: the floor bounds the degenerate species, and only them", {
  O <- log(trichoptera$Offset)
  X <- model.matrix(~ 1, trichoptera)
  f <- Abundance ~ 1 + offset(log(Offset))
  floor <- 1e-1 # high enough to bind on these data

  ## no degenerate species on these data: the floor does nothing
  free  <- PLNnetwork(f, trichoptera, control = PLNnetwork_param(trace = 0, n_penalties = 5))
  quiet <- PLNnetwork(f, trichoptera, control = PLNnetwork_param(trace = 0, n_penalties = 5, latent_floor = floor))
  expect_equal(vapply(quiet$models, function(m) m$optim_par$n_floor, numeric(1)), rep(0, 5))
  expect_true(all(lengths(lapply(quiet$models, function(m) m$floored_species)) == 0))
  expect_equal(quiet$criteria$n_edges, free$criteria$n_edges)
  expect_equal(quiet$criteria$loglik, free$criteria$loglik, tolerance = 1e-3)
  expect_length(free$models[[1]]$floored_species, 0)

  ## a low threshold stands for degenerate species
  old <- options(PLNmodels.latent_variance_threshold = 3)
  on.exit(options(old))
  nets <- suppressWarnings( # the degenerate species are reported in a warning
    PLNnetwork(f, trichoptera, control = PLNnetwork_param(trace = 0, n_penalties = 5, latent_floor = floor))
  )
  floored <- lapply(nets$models, function(m) m$floored_species)
  expect_gt(length(floored[[5]]), 0)
  expect_lt(length(floored[[5]]), nets$models[[5]]$p)
  ## a species bounded at a penalty stays so down the path
  for (k in 2:5) expect_true(all(floored[[k - 1]] %in% floored[[k]]))
  for (m in nets$models) {
    expect_equal(m$latent_floor, floor)
    counts <- exp(O + m$var_par$M)
    bounded <- colnames(trichoptera$Abundance) %in% m$floored_species
    ## the bounded species are above the floor, the others are free to go below
    if (any(bounded)) expect_gte(min(counts[, bounded]), floor * (1 - 1e-10))
    expect_equal(m$optim_par$n_floor, sum(counts[, bounded] <= floor * (1 + 1e-10)))
    ## the lower bound is the one of the parameters returned
    expect_equal(m$loglik, sum(PLNmodels:::elbo_fixed_precision(
      as.matrix(trichoptera$Abundance), X, matrix(O, nrow(X), m$p), m$model_par$B,
      m$var_par$M, m$var_par$S2, m$model_par$Omega)), tolerance = 1e-8)
  }
  last <- nets$models[[5]]
  expect_gt(last$optim_par$n_floor, 0)
  expect_lt(min(exp(O + last$var_par$M)[, !colnames(trichoptera$Abundance) %in% last$floored_species]), floor)
})

test_that("the projection on the floor sets the variance of the bounded cells to its optimum", {
  set.seed(2)
  n <- 12; p <- 4
  data <- list(Y = matrix(rpois(n * p, 2), n, p), X = matrix(1, n, 1), O = matrix(0, n, p), w = rep(1, n))
  Omega <- diag(runif(p, 0.5, 2))
  out <- list(M = matrix(rnorm(n * p, -3, 2), n, p), S2 = matrix(0.5, n, p), B = matrix(0, 1, p))
  floor <- 0.05
  below <- out$M < log(floor)
  proj <- PLNmodels:::project_latent_floor(out, data, Omega, floor)
  expect_equal(proj$n_floor, sum(below))
  expect_true(all(proj$M >= log(floor)))
  expect_equal(proj$M[!below], out$M[!below])
  expect_equal(proj$S2[!below], out$S2[!below])
  ## 1 / S2 = Omega_jj + exp(O + M + S2 / 2) on the bounded cells
  om <- matrix(diag(Omega), n, p, byrow = TRUE)
  expect_equal(1 / proj$S2[below], (om + exp(proj$M + proj$S2 / 2))[below], tolerance = 1e-8)
  ## B, Sigma and the bound are those of the projected parameters
  expect_equal(as.vector(proj$B), colMeans(proj$M))
  R <- sweep(proj$M, 2, colMeans(proj$M))
  expect_equal(proj$Sigma, crossprod(R) / n + diag(colMeans(proj$S2)))
  expect_equal(proj$Ji, PLNmodels:::elbo_fixed_precision(data$Y, data$X, data$O, proj$B, proj$M, proj$S2, Omega))

  ## restricted to some species, the others are left as they are
  some <- c(TRUE, FALSE, TRUE, FALSE)
  part <- PLNmodels:::project_latent_floor(out, data, Omega, floor, some)
  expect_equal(part$M[, !some], out$M[, !some])
  expect_equal(part$S2[, !some], out$S2[, !some])
  expect_true(all(part$M[, some] >= log(floor)))
  expect_equal(part$n_floor, sum(below[, some]))
})

test_that("ZIPLNnetwork and stability selection follow the scale of the penalty", {
  tri <- trichoptera[1:25, ]
  tri$Abundance <- tri$Abundance[, 1:6]
  zi <- ZIPLNnetwork(Abundance ~ 1, tri, control = ZIPLNnetwork_param(trace = 0, n_penalties = 4, penalty_scale = "correlation"))
  expect_true(all(zi$penalties > 0 & zi$penalties <= 1))
  expect_equal(zi$models[[1]]$penalty_scale, "correlation")
  expect_true(all(is.finite(zi$criteria$loglik)))

  nets <- PLNnetwork(Abundance ~ 1, tri,
                     control = PLNnetwork_param(trace = 0, n_penalties = 4, penalty_scale = "correlation", latent_floor = 1e-3))
  subs <- replicate(2, sample.int(nrow(tri), 15), simplify = FALSE)
  capture.output(nets$stability_selection(subsamples = subs))
  expect_false(is.null(nets$stability_path))
})

test_that("a species absent from a group of samples is not a hub on the correlation scale", {
  skip_on_cran()
  ## a sparse network on 20 species, 60 samples; two species are then made
  ## absent from 65 % of the samples, and the model is fitted without the
  ## covariate that would explain it
  set.seed(3)
  n <- 60; p <- 20
  A <- matrix(0, p, p); A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, 0.1); A <- A + t(A)
  Omega <- 0.5 * A; diag(Omega) <- rowSums(abs(Omega)) + 0.1
  Z <- matrix(rnorm(n * p), n, p) %*% chol(cov2cor(solve(Omega))) + matrix(runif(p, 2, 4), n, p, byrow = TRUE)
  Y <- matrix(rpois(n * p, exp(Z)), n, p, dimnames = list(paste0("s", 1:n), paste0("sp", 1:p)))
  contaminated <- c(3, 11)
  for (j in contaminated) Y[sample(n, round(0.65 * n)), j] <- 0
  dat <- suppressMessages(prepare_data(Y, data.frame(id = 1:n, row.names = rownames(Y)), offset = "none"))

  ## share of the edges touching a contaminated species, in the model with the true number of edges
  contaminated_share <- function(nets) {
    n_edges <- vapply(nets$models, function(m) m$n_edges, numeric(1))
    net <- as.matrix(nets$models[[which.min(abs(n_edges - sum(A) / 2))]]$latent_network("support"))
    (sum(net) / 2 - sum(net[-contaminated, -contaminated]) / 2) / (sum(net) / 2)
  }
  covariance  <- suppressWarnings(PLNnetwork(Abundance ~ 1, dat, control = PLNnetwork_param(trace = 0, min_ratio = 0.05)))
  correlation <- suppressWarnings(PLNnetwork(Abundance ~ 1, dat, control = PLNnetwork_param(trace = 0, min_ratio = 0.05, penalty_scale = "correlation")))
  ## 2 species out of 20: about 20 % of the edges would touch them by chance
  expect_gt(contaminated_share(covariance), 0.5)
  expect_lt(contaminated_share(correlation), 0.3)
})
