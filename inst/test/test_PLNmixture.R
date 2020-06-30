library(PLNmodels)
library(aricode)

data("iris")
n <- nrow(iris)
p <- 4
lambda <- exp(as.matrix(iris[, 1:p]))
count <- matrix(rpois (n * p, lambda), n, p)
myPLN <- PLNmixture(count ~ 1, clusters = 1:5, control_main = list(covariance = "spherical"))
myPLN$plot()
myPLN$plot_objective()
# bestModel <- myPLN$getBestModel("ICL")
bestModel <- myPLN$getModel(3)

aricode::ARI(bestModel$memberships, iris$Species)
pairs(iris[, 1:4], col = bestModel$memberships, pch = as.numeric(iris$Species) )

kmeans_cl <- kmeans(log(1 + count), centers = 3)$cl
# kmeans_cl <- kmeans(count, centers = 3)$cl
aricode::ARI(kmeans_cl, iris$Species)
pairs(iris[, 1:4], col = kmeans_cl, pch = as.numeric(iris$Species) )

aricode::ARI(bestModel$memberships, kmeans_cl)

data("mollusk")
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)
mollusk_mixture <- PLNmixture(Abundance ~ 1 +  offset(log(Offset)), clusters = 1:6, data = mollusc, control_main = list(covariance = "spherical", cores = 3))
plot(mollusk_mixture$models[[5]]$optim_par$objective, type = "l")

mollusk_mixture$plot()
bestModel <- mollusk_mixture$getModel(4)
plot(bestModel$optim_par$convergence, type = "l")

mollusk_PCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), ranks = 1:5, data = mollusc, control_main = list(cores = 5))
plot(getBestModel(mollusk_PCA), map = "individual", ind_cols = mollusc$site)
plot(getBestModel(mollusk_PCA), map = "individual", axes = c(1,2), ind_cols = as.factor(mollusk_mixture$models[[4]]$memberships))

