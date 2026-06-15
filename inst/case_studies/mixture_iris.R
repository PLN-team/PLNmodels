data(iris)

library(PLNmodels)
library(tidyverse)
library(viridisLite)

count <- iris %>% dplyr::select(-Species) %>% exp() %>% round() %>% as.matrix()
covariate <- data.frame(Species = iris$Species)
iris_data <- prepare_data(count, covariate)
my_mixtures <-  PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:5, data = iris_data)

plot(my_mixtures)

myPLN <- getBestModel(my_mixtures)
myPLN <- getModel(my_mixtures, 3)
plot(myPLN, type = "pca")
plot(myPLN, type = "matrix")

aricode::ARI(myPLN$memberships, iris$Species)

my_mixtures_covar <-  PLNmixture(Abundance ~ 0 + Species + offset(log(Offset)), clusters = 1:3, data = iris_data)

plot(my_mixtures_covar)

myPLN_covar <- getModel(my_mixtures_covar, 2)
plot(myPLN_covar, "pca")
plot(myPLN_covar, "matrix")
