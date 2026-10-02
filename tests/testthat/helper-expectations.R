## Expectations shared by several test files

## A collection of network fits (PLNnetwork or ZIPLNnetwork) of class `fit_class`:
## S3 methods, fields and R6 methods; `...` goes to the stability selection
expect_network_family <- function(family, fit_class, ...) {
  ## S3 methods
  expect_true(PLNmodels:::isNetworkfamily(family))
  expect_is(plot(family), "ggplot")
  expect_is(plot(family, reverse = TRUE), "ggplot")
  expect_is(plot(family, type = "diagnostic"), "ggplot")
  expect_is(getBestModel(family), fit_class)
  expect_is(getModel(family, family$penalties[1]), fit_class)

  ## Field access
  expect_true(all(family$penalties > 0))
  expect_null(family$stability_path)
  expect_true(anyNA(family$stability))

  ## Other R6 methods
  expect_true(is.data.frame(family$coefficient_path()))
  n <- nrow(family$responses)
  subs <- replicate(2, sample.int(n, size = n/2), simplify = FALSE)
  family$stability_selection(subsamples = subs, ...)
  expect_is(plot(family, type = "stability"), "ggplot")
  expect_true(!is.null(family$stability_path))
  expect_true(inherits(family$plot(), "ggplot"))
  expect_true(inherits(family$plot_objective(), "ggplot"))
  expect_true(inherits(family$plot_stars(), "ggplot"))
}

## Univariate fits (lists, one fit per variable, for each covariance model) are
## equivalent to each other and to the multivariate diagonal fit
expect_univariate_equivalence <- function(univariate_full, univariate_diagonal, univariate_spherical,
                                          multivariate_diagonal) {
  univariate <- list(univariate_full, univariate_diagonal, univariate_spherical)
  variances  <- purrr::map(univariate, ~ purrr::map_dbl(purrr::map(.x, sigma), as.double))
  for (fits in univariate) {
    expect_true(all.equal(purrr::map_dbl(fits, "nb_param"), purrr::map_dbl(univariate_full, "nb_param")))
    expect_true(all.equal(sum(purrr::map_dbl(fits, "loglik")), multivariate_diagonal$loglik, tolerance = 1e-2))
  }
  expect_true(all.equal(variances[[3]], variances[[2]], tolerance = .25))
  expect_true(all.equal(variances[[3]], variances[[1]], tolerance = .25))
  expect_true(all.equal(variances[[2]], variances[[1]], tolerance = .25))
}
