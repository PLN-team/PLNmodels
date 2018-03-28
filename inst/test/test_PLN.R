library(PLNmodels)
library(ggplot2)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau) ## must be a matrix

profiling_1 <- profr::profr(model_1 <- PLN(abundance ~ 1, control = list(newpar = FALSE)))
profiling_2 <- profr::profr(model_2 <- PLN(abundance ~ 1, control = list(newpar = TRUE, method = "CCSAQ")))


library(microbenchmark)
res <- microbenchmark(nloptr = PLN(abundance ~ 1, control = list(trace = FALSE, newpar = FALSE)),
                      nlopt  = PLN(abundance ~ 1, control = list(trace = FALSE, newpar = TRUE , method = "CCSAQ")), times = 20)

par(mfrow = c(1,3))
plot(profiling_1)
plot(profiling_2)
plot(res, log = "y", las = 3)
