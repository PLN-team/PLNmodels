library(PLNmodels)
library(ade4)
data("trichometeo")
abundance <- as.matrix(trichometeo$fau) ## must be a matrix

model_1 <- PLN(abundance ~ 1, control = list(newpar = FALSE))
model_2 <- PLN(abundance ~ 1, control = list(newpar = TRUE ))

library(microbenchmark)
res <- microbenchmark(nloptr = PLN(abundance ~ 1, control = list(trace = FALSE, newpar = FALSE)),
                      nlopt  = PLN(abundance ~ 1, control = list(trace = FALSE, newpar = TRUE )), times = 20)
library(ggplot2)
autoplot(res)
