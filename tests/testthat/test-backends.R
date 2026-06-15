context("test-backends")

## Non-regression guards for backend correctness.
## Key reference values come from benchmark runs on code-enhancement branch.

data(barents)

## --- PLN barents : guard against XTOL premature convergence bug (master had ll ~ -8520) ---
test_that("PLN builtin on barents converges well above the XTOL-bug threshold", {
  m <- PLN(Abundance ~ Depth + Temperature + offset(log(Offset)),
           data = barents,
           control = PLN_param(backend = "builtin", trace = 0))
  ## master (nlopt) was stuck at -8520 due to X-scaling / XTOL bug.
  ## CE builtin reaches ~ -4402; we guard with a conservative threshold of -6000.
  expect_gt(m$loglik, -6000)
})

test_that("PLN nlopt on barents also surpasses the XTOL-bug threshold", {
  m <- PLN(Abundance ~ Depth + Temperature + offset(log(Offset)),
           data = barents,
           control = PLN_param(backend = "nlopt", trace = 0))
  expect_gt(m$loglik, -6000)
})

## --- builtin reaches at least as good loglik as nlopt on barents ---
test_that("PLN builtin loglik >= nlopt loglik on barents (full covariance)", {
  ctrl_b <- PLN_param(backend = "builtin", trace = 0)
  ctrl_n <- PLN_param(backend = "nlopt",   trace = 0)
  ll_b <- PLN(Abundance ~ Depth + Temperature + offset(log(Offset)), data = barents, control = ctrl_b)$loglik
  ll_n <- PLN(Abundance ~ Depth + Temperature + offset(log(Offset)), data = barents, control = ctrl_n)$loglik
  ## builtin should find an equal or better solution; allow 1 unit tolerance for stochasticity
  expect_gte(ll_b, ll_n - 1)
})
