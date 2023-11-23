
## data
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
nb_blocks <- 1:15

## Start with PLN to get back latent position
myPLN <- PLN(Abundance ~ 1,  data = trichoptera)

## Kmean init on the variaitonal means
Means <- t(myPLN$var_par$M)
init_kmeans_means <- lapply(nb_blocks, function(k) {
  kmeans(Means, centers = k, nstart = 30)$cl
})
blocks1 <- PLNblock(Abundance ~ 1,  data = trichoptera, nb_blocks = 1:15,
                    control = PLNblock_param(inception = myPLN, init_cl = init_kmeans_means))

D <- 1 - cov2cor(myPLN$model_par$Sigma)
init_ward2_cor <- hclust(as.dist(D), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
blocks2 <- PLNblock(Abundance ~ 1,  data = trichoptera, nb_blocks = 1:15,
                    control = PLNblock_param(inception = myPLN, init_cl = init_ward2_cor))

Sbar <- colSums(myPLN$var_par$S2)
D <- sqrt(as.matrix(dist(t(myPLN$var_par$M)^2)) + outer(Sbar,rep(1,myPLN$p)) + outer(rep(1, myPLN$p), Sbar))
init_kmeans_dist <- lapply(nb_blocks, function(k) {
  kmeans(D, centers = k, nstart = 30)$cl
})
blocks3 <- PLNblock(Abundance ~ 1,  data = trichoptera, nb_blocks = 1:15,
                    control = PLNblock_param(inception = myPLN, init_cl = init_kmeans_dist))

init_ward2_dist <- hclust(as.dist(D), method = "ward.D2") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
blocks4 <- PLNblock(Abundance ~ 1,  data = trichoptera, nb_blocks = 1:15,
                    control = PLNblock_param(inception = myPLN, init_cl = init_ward2_dist))
