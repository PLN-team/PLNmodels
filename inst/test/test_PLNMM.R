library(PLNmodels)
library(aricode)

data("iris")

lambda <- exp(as.matrix(iris[, 1:4]))
count <- matrix(rpois (nrow(iris) * 4, lambda), nrow(iris), 4)

myPLN <- PLNMM(count ~ 1, clusters = 1:4)
myPLN$plot()
bestModel <- myPLN$getBestModel("ICL")
plot(bestModel$optim_par$convergence, type = "l")
aricode::ARI(bestModel$memberships, iris$Species)
pairs(iris[, 1:4], col = bestModel$memberships, pch = as.numeric(iris$Species) )

kmeans_cl <- kmeans(log(1 + count), centers = 3)$cl
aricode::ARI(kmeans_cl, iris$Species)
pairs(iris[, 1:4], col = kmeans_cl, pch = as.numeric(iris$Species) )
