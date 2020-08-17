# une vue synthétique (genre heatmap) de la matrice de comptage prédite
# (A, mais sans l'offset et sans les covariables) avec les échantillons
# triés par groupe, pour faire ressortir les blocs d'échantillons et
# espérer qu'il y ait des blocs de taxa associés (mais c'est mon biais LBM
# qui parle)

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

my_mixtures <-  PLNmixture(Abundance ~ 1 + Species + offset(log(Offset)), clusters = 1:3, data = iris_data, control_main = list(core = nb_cores))

