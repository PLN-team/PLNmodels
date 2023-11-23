library(tidyverse)
library(viridis)

## setting up future for parallelism
nb_cores <- 10
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

## data
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
nb_blocks <- 1:15

data(oaks)
nb_blocks <- 3:70

data <- oaks

## Start with PLN to get back latent position
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)),  data = data)

## Kmean init on the variaitonal means
Means <- t(myPLN$var_par$M)
D <- as.matrix(dist(t(myPLN$var_par$M)^2))

init_km_1 <- lapply(nb_blocks, function(k) {
  kmeans(Means, centers = k, iter = 100, nstart = 30)$cl
})
blocks_km_1 <- PLNblock(Abundance ~ 1 + offset(log(Offset)),data = data, nb_blocks = nb_blocks,
                    control = PLNblock_param(inception = myPLN, init_cl = init_km_1))

init_km_2 <- lapply(nb_blocks, function(k) {
  kmeans(D, centers = k, iter = 100,  nstart = 30)$cl
})
blocks_km_2 <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                    control = PLNblock_param(inception = myPLN, init_cl = init_km_2))

init_ward2_1 <- hclust(as.dist(D), method = "ward.D2") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
blocks_hc_1 <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                    control = PLNblock_param(inception = myPLN, init_cl = init_ward2_1))

D <- 1 - cov2cor(myPLN$model_par$Sigma)
init_ward2_2 <- hclust(as.dist(D), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
blocks_hc_2 <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                    control = PLNblock_param(inception = myPLN, init_cl = init_ward2_2))

init_ward2_3 <- hclustvar(myPLN$var_par$M) %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
blocks_hc_3 <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                        control = PLNblock_param(inception = myPLN, init_cl = init_ward2_3))

bind_rows(bind_cols(blocks_km_1$criteria, method = "kmeans on M"),
          bind_cols(blocks_km_2$criteria, method = "kmeans on Euclidian distances on M "),
          bind_cols(blocks_hc_1$criteria, method = "Ward D2 on Euclidian distances"),
          bind_cols(blocks_hc_2$criteria, method = "Ward D2 on correlation"),
          bind_cols(blocks_hc_3$criteria, method = "ClustOfVar")) %>%
  group_by(method) %>%
  ggplot() + aes(x = param, y = loglik, colour = method) + geom_line() + scale_color_viridis_d() + theme_bw()

future::plan("sequential")
