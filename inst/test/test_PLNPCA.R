library(PLNmodels)
library(ggplot2)
library(profr)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau)

profiling <- profr::profr(model <- PLNPCA(abundance ~ 1, ranks = 1:5))
plot(profiling)

library(microbenchmark)
res <- microbenchmark(nloptr = PLNPCA(abundance ~ 1, ranks = 1:5), times = 20)
autoplot(res)
