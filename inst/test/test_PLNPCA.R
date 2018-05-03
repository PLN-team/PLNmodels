library(PLNmodels)
library(ggplot2)
library(testthat)
library(profr)
library(ade4)
data("trichometeo")

abundance <- as.matrix(trichometeo$fau)

profiling1 <- profr::profr(model1 <- PLNPCA(abundance, ranks = 1:5))
profiling2 <- profr::profr(model2 <- PLNPCA(abundance ~ 1, ranks = 1:5))

expect_equivalent(model1, model2)

par(mfrow = c(2,1))
plot(profiling1)
plot(profiling2)

library(microbenchmark)
res <- microbenchmark(PLNPCA = PLNPCA(abundance, ranks = 1:5, control.main = list(trace = 0)), times = 10)
autoplot(res)
