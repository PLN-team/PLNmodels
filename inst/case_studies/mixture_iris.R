data(iris)

library(PLNmodels)
library(tidyverse)
library(viridisLite)

nb_cores <- 10

## Dimension reduction with PCA

count <- iris %>% dplyr::select(-Species) %>% exp() %>% round() %>% as.matrix()
covariate <- data.frame(Species = iris$Species)
iris_data <- prepare_data(count, covariate)
my_mixtures <-  PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:5, data = iris_data, control_main = list(core = nb_cores))

plot(my_mixtures)

myPLN <- getBestModel(my_mixtures)
myPLN$plot_clustering_pca()
myPLN$plot_clustering_data()

aricode::ARI(myPLN$memberships, iris$Species)

my_mixtures_covar <-  PLNmixture(Abundance ~ 1 + Species + offset(log(Offset)), clusters = 1:3, data = iris_data, control_main = list(core = nb_cores))

plot(my_mixtures_covar)

myPLN_covar <- getBestModel(my_mixtures_covar)
myPLN_covar$plot_clustering_pca()
myPLN_covar$plot_clustering_data()
