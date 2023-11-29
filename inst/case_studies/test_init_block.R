library(PLNmodels)
library(ClustOfVar)
library(mclust)
library(tidyverse)
library(viridis)

## setting up future for parallelism
nb_cores <- 20
options(future.fork.enable = TRUE)
future::plan("multicore", workers = nb_cores)

## data
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
nb_blocks <- 2:15
data <- trichoptera

data(oaks)
nb_blocks <- seq(3, 72, by = 3)
data <- oaks

elbo <- function(clusters, PLN_fit, data, formula) {
  Tau <- PLNmodels:::as_indicator(clusters)
  mf <- model.frame(formula, data)
  Y <- model.response(mf)
  X <- model.matrix(mf, data)
  O <- model.offset(mf)
  n <- nrow(Y)

  B  <- PLN_fit$model_par$B
  M  <- PLN_fit$var_par$M %*% Tau
  S2 <- (PLN_fit$var_par$S %*% Tau)**2
  Sigma <- (1/n) * (crossprod(M) + diag(colSums(S2)))

  A <- exp(O + X %*% B) * exp(M + .5 * S2) %*% t(Tau);
  objective <- sum((A - Y * (O + X %*% B + M %*% t(Tau)))) - 0.5 * sum(log(S2)) + 0.5 * n * determinant(Sigma)$modulus
  -objective
}

## Start with PLN to get back latent position
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = data)

formula <- Abundance ~ 1 + offset(log(Offset))
## Kmean init on the variatonal means
Mt <- t(myPLN$var_par$M)
n <- ncol(Mt)
init_boot <- future_replicate(100, {
  Mts <- Mt[, sample.int(n)]
  sapply(nb_blocks, function(k) {kmeans(Mts, centers = k, iter = 100, nstart = 30)$cl})
}, simplify = FALSE)
res <- sapply(init_boot, function(samp) apply(samp, 2, elbo, myPLN, data, formula))
idx <- apply(res, 1, which.max)
init_km_boot <- vector("list", length(idx))
for (i in seq_along(idx)){
  init_km_boot[[i]] <- init_boot[[idx[i]]][, i]
}

init_km <- lapply(nb_blocks, function(k) {kmeans(Mt, centers = k, iter = 100, nstart = 30)$cl})
init_hc <- hclust(as.dist(1 - cov2cor(myPLN$model_par$Sigma)), method = "complete") %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
init_hcvar <- hclustvar(myPLN$var_par$M) %>% cutree(nb_blocks) %>% as.data.frame() %>% as.list()
init_gmm <- lapply(nb_blocks, function(k) Mclust(Mt, G = k, modelNames = c( "EEI"), verbose = FALSE)$classification)

idx <- apply(data.frame(
  crit_km = sapply(init_km, elbo, myPLN, data, formula),
  crit_hc =  sapply(init_hc, elbo, myPLN, data, formula),
  crit_hcvar = sapply(init_hcvar, elbo, myPLN, data, formula),
  crit_gmm   = sapply(init_gmm, elbo, myPLN, data, formula)), 1, which.max
)
all_init <- list(init_km, init_hc, init_hcvar, init_gmm)[idx]
init_cs <- vector("list", length(idx))
for (i in seq_along(idx)){
  init_cs[[i]] <- all_init[[idx[i]]][[i]]
}

x11();
data.frame( blocks = rep(nb_blocks, 4),
            criteria = c(crit_km, crit_hc, crit_hcvar, crit_gmm),
            method = rep(c("kmeans", "complete linkage", "ClustOfVar", "Spherical GMM"), each = length(nb_blocks))) %>%
  ggplot() + aes(x = blocks, y = criteria, colour = method, group = method) +
  geom_line() + coord_trans(y="log10") + scale_color_viridis_d() + theme_bw()

blocks_gmm <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                       control = PLNblock_param(inception = myPLN, init_cl = init_gmm))
blocks_hcvar <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                         control = PLNblock_param(inception = myPLN, init_cl = init_hcvar))
blocks_hc <- PLNblock(Abundance ~ 1 + offset(log(Offset)), data = data, nb_blocks = nb_blocks,
                      control = PLNblock_param(inception = myPLN, init_cl = init_hc))
blocks_km <- PLNblock(Abundance ~ 1 + offset(log(Offset)),data = data, nb_blocks = nb_blocks,
                      control = PLNblock_param(inception = myPLN, init_cl = init_km))
blocks_cs <- PLNblock(Abundance ~ 1 + offset(log(Offset)),data = data, nb_blocks = nb_blocks,
                      control = PLNblock_param(inception = myPLN, init_cl = init_cs))
blocks_km_boot <- PLNblock(Abundance ~ 1 + offset(log(Offset)),data = data, nb_blocks = nb_blocks,
                      control = PLNblock_param(inception = myPLN, init_cl = init_km_boot))

x11()
bind_rows(bind_cols(blocks_km$criteria, method = "kmeans"),
          bind_cols(blocks_km_boot$criteria, method = "kmeans bootstrap"),
          bind_cols(blocks_hc$criteria, method = "HC (complete linkage)"),
          bind_cols(blocks_hcvar$criteria, method = "ClustOfVar"),
          bind_cols(blocks_gmm$criteria, method = "GMM"),
          bind_cols(blocks_cs$criteria, method = "Consensus")) %>%
  group_by(method) %>%
  ggplot() + aes(x = param, y = loglik, colour = method) + geom_line() + scale_color_viridis_d() + theme_bw()

future::plan("sequential")
