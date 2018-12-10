library(PLNmodels)
library(testthat)

## timings (full, spherical, diagonal)
res <- microbenchmark::microbenchmark(
  spherical = PLN(Abundance ~ 1, data = trichoptera, control = list(covariance  ="spherical", trace = 0)),
  diagonal  = PLN(Abundance ~ 1, data = trichoptera, control = list(covariance  ="diagonal", trace = 0)),
  full      = PLN(Abundance ~ 1, data = trichoptera, control = list(trace = 0)),
  times = 20
)
summary(res)

# res <- microbenchmark::microbenchmark(
#   uw = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical",  control = list(trace = 0)),
#   ## equivalent weigths
#   w  = PLN(Abundance ~ 1, data = trichoptera, covariance  ="spherical", weights = rep(1.0, nrow(trichoptera)), control = list(trace = 0)),
#   times = 20
# )
# summary(res)
