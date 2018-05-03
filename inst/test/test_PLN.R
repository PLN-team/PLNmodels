library(PLNmodels)
library(ggplot2)
library(testthat)
library(profr)
library(ade4)
data("trichometeo")

abundance <- as.matrix(trichometeo$fau)

profiling1 <- profr::profr(model1 <- PLN(abundance))
profiling2 <- profr::profr(model2 <- PLN(abundance ~ 1))

expect_equivalent(model1, model2)

par(mfrow = c(2,1))
plot(profiling1)
plot(profiling2)

library(microbenchmark)
res <- microbenchmark(nloptr = PLN(abundance, control = list(trace = 0)), times = 20)
autoplot(res)
