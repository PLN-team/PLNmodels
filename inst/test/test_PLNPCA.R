library(PLNmodels)
library(ggplot2)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau)

profiling_1 <- profr::profr(model_1 <- PLNPCA(abundance ~ 1, ranks = 1:5, control.init = list(nloptr = TRUE), control.main = list(nloptr = TRUE )))
profiling_2 <- profr::profr(model_2 <- PLNPCA(abundance ~ 1, ranks = 1:5, control.main = list(nloptr = FALSE)))

library(microbenchmark)
res <- microbenchmark(nloptr = PLNPCA(abundance ~ 1, ranks = 1:5, control.init = list(nloptr = TRUE), control.main = list(trace = FALSE, nloptr = TRUE )),
                      nlopt  = PLNPCA(abundance ~ 1, ranks = 1:5, control.main = list(trace  = FALSE)), times = 20)

par(mfrow = c(1,3))
plot(profiling_1)
plot(profiling_2)
plot(res, log = "y", las = 3)
