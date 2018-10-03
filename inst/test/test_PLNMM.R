library(PLNmodels)
library(aricode)

data("iris")

lambda <- exp(as.matrix(iris[, 1:4]))
count <- matrix(rpois (nrow(iris) * 4, lambda), nrow(iris), 4)

myPLN <- PLNMM(count ~ 1, clusters = 2:4)

bestModel <- myPLN$getBestModel('BIC')

aricode::ARI(bestModel$memberships, iris$Species)

pairs(iris[, 1:4], col = bestModel$memberships, pch = as.numeric(iris$Species) )

plot(bestModel$optim_par$convergence, type = "l")
