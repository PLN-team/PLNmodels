library(PLNmodels)
library(ggplot2)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau) ## must be a matrix

profiling <- profr::profr(model <- PLN(abundance ~ 1))
plot(profiling)

library(microbenchmark)
res <- microbenchmark(nloptr = PLN(abundance ~ 1), times = 20)
autoplot(res)
