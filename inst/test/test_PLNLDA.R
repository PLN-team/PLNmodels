library(PLNmodels)
library(ggplot2)
library(testthat)
library(microbenchmark)
library(profr)
library(ade4)
data("trichometeo")

abundance <- as.matrix(trichometeo$fau)
night_grp <- as.character(trichometeo$cla)

profiling1 <- profr::profr(model1 <- PLNLDA(abundance, night_grp))
profiling2 <- profr::profr(model2 <- PLNLDA(abundance ~ 0, night_grp))

expect_equivalent(model1, model2)

par(mfrow = c(2,1))
plot(profiling1)
plot(profiling2)

res <- microbenchmark(PLNLDA = PLNLDA(abundance, night_grp, control = list(trace = 0)), times = 20)
autoplot(res)
