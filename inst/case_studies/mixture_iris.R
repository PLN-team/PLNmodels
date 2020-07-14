# 3/ l'inévitable scatterplot des positions dans l'espace latent
# (éventuellement après avoir fait une ACP sur M) avec coloration en
# fonction du groupe du MAP
# 4/ une vue synthétique (genre heatmap) de la matrice de comptage prédite
# (A, mais sans l'offset et sans les covariables) avec les échantillons
# triés par groupe, pour faire ressortir les blocs d'échantillons et
# espérer qu'il y ait des blocs de taxa associés (mais c'est mon biais LBM
# qui parle)

data(iris)

library(PLNmodels)
library(tidyverse)
library(viridisLite)
library(viridisLite)

nb_cores <- 4

## Dimension reduction with PCA

count <- iris %>% select(-Species) %>% exp() %>% round() %>% as.matrix()
covariate <- data.frame(iris$Species)
iris_data <- prepare_data(count, covariate)
my_mixtures <-  PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:5, data = iris_data, control_main = list(covariance = 'spherical', core = nb_cores))

plot(my_mixtures)

myPLN <- my_mixtures %>% getBestModel()
myPLN$plot_clustering_pca()

aricode::ARI(myPLN$memberships, iris$Species)

p <- Reduce("+", Map(function(pi, comp) {pi * comp$var_par$M}, myPLN$mixtureParam, myPLN$components)) %>%
  as_tibble() %>% setNames(c("sepal_length", "sepal_width", "petal_length", "petal_width")) %>%
  add_column(memberships = factor(myPLN$memberships)) %>%
  GGally::ggpairs(aes(colour = memberships), upper = list(continuous = "density")) + theme_bw()
p
