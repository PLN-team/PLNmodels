library(PLNmodels)
library(parallel)
library(purrr)

data(oaks)
all_mixtures <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = oaks, clusters = 1:7, control_main = list(iterates = 0))

k <- 3

my_model <- all_mixtures$models[[k]]

my_model$entropy
my_model$entropy_clustering
my_model$entropy_latent

k <- 7
tau_candidates <- lapply(combn(k, 2, simplify = FALSE), function(couple) {
  tau <- my_mixtures$models[[k]]$posteriorProb
  i <- min(couple)
  j <- max(couple)
  tau_merged <- tau[, -j]
  tau_merged[, i] <- rowMeans(tau[, c(i,j)])
  tau_merged
})

control <- PLNmodels:::PLNmixture_param(list(maxit_out = 2), my_model$n, my_model$p)

loglik_candidates <- mclapply(tau_candidates, function(tau) {

  model <- PLNmodels:::PLNmixturefit$new(all_mixtures$responses, all_mixtures$covariates, all_mixtures$offsets, tau, as.formula(Abundance ~ 1 + offset(log(Offset))), NULL, control)
  model$optimize(all_mixtures$responses, all_mixtures$covariates, all_mixtures$offsets, control)

  model$loglik
})

