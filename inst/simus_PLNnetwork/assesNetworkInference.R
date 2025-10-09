working.dir <- "~/Desktop/PLNnetwork-simu/"

library(tidyverse)
library(viridis)
library(pbmcapply)
source(paste0(working.dir,"func_data_generation.R"))
source(paste0(working.dir,"func_networks.R"))
source(paste0(working.dir,"func_utils.R"))

setwd(working.dir)

# methods considered
methods_names <- set_names(
  c("SPIEC-EASI",
    "sparse PLN",
    # "sparse PLNCov",
    "sparse ZIPLNCov",
    "sparse ZIPLNCol",
    "GLasso",
    "NeighSelec",
    "SparseCC"
    ),
  nm = c("spiecEasi_network",
         "sparsePLN_network",
         # "sparsePLNCov_network",
         "sparseZIPLNCov_network",
         "sparseZIPLNCol_network",
         "graphical_lasso_network",
         "neighborhood_selection_network",
         "sparCC_network"
    )
  )

## Simulation parameters that remain fixed
p <- 50
k <- 25
offset_effect <- 2  ## the larger, the less effect
thres_auc <- .25

seq_n  <- c(25, 50, 100)
seq_covar <- 1
seq_covar_zi <- "contrasted" #"skew-right" # "skew-left" "contrasted"
cores <- 5

params <- expand.grid(seq_n, seq_covar, seq_covar_zi)
colnames(params) <- c("n", "covar_effect", "covar_zi")
nsimu <- 50

one_simu <- function(n, covar_effect, zi_effect, network) {
  Sigma_star <- graph2cov(network, v = 0.3, u = 0.1)
  group <- factor(sample(1:3, n, replace = TRUE))
  X <- cbind(1, rnorm(n, 2, 1) #  model.matrix(~ group + 0)
  mu <- X %*% matrix(runif(ncol(X)*p, -covar_effect,  covar_effect) , ncol(X), p)

  par_zi <- switch(zi_effect, "contrasted" = c(-5,5), "skew-right" = c(0,5), "skew-left" = c(-5,0))
  group0 <- factor(sample(1:2, n, replace = TRUE))
  X0 <- model.matrix(~ group0 + 0)
  mu0 <- X0 %*% matrix(runif(ncol(X0)*p, par_zi[1],  par_zi[2]) , ncol(X0), p)
  
  N  <- rnegbin(n, mu = 1000, theta = offset_effect); N[N == 0] <- 1000
  Y <- rMLN(n, mu = mu, Sigma = Sigma_star, N = N)
  
  W <- matrix(rbinom(n*p, size = 1, prob = sigmoid(mu0)), n, p)
  Y <- Y * W
  which_zero  <- which(colSums(Y) == 0)
  for (j in which_zero) Y[sample.int(n, 1), j] <- 1

  # run inferences
  methods_call    <- lapply(names(methods_names), function(method) call(method, Y , X, X0))
  methods_results <- lapply(methods_call, eval)

  # performance
  res_roc <- lapply(methods_results, perf_roc, as_adj(network))
  names(res_roc) <- methods_names
  res_auc <- sapply(res_roc, perf_auc, thres_auc)
  names(res_auc) <- methods_names
  res_aupr <- sapply(res_roc, perf_aupr)
  names(res_aupr) <- methods_names
  
  tidy_res <- data.frame(
    fallout = unlist(lapply(res_roc, function(res) res$fallout )),
    recall  = unlist(lapply(res_roc, function(res) res$recall )) ,
    precision  = unlist(lapply(res_roc, function(res) res$precision )) ,
    pAUC      = rep(round(res_auc, 3), sapply(res_roc, function(res) length(res$fallout))),
    AUPR      = rep(round(res_aupr, 3), sapply(res_roc, function(res) length(res$fallout))),
    method = rep(methods_names, sapply(res_roc, function(res) length(res$fallout))),
    method_auc = rep(paste(methods_names, round(res_auc, 3), sep = " - pAUC = "), sapply(res_roc, function(res) length(res$fallout))),
    row.names = NULL
  )
  tidy_res$sample_size  <- switch(as.character(n), "25" = "n = p/2", "50" = "n = p", "100" = "n = 2p")
  tidy_res$covar_effect <- switch(as.character(covar_effect), "1" = "small", "2" = "medium", "3" = "large")
  tidy_res$covar_zi     <- zi_effect
  tidy_res
}

# out <- one_simu(100, 1, "contrasted", network = sample_er(p, k))

res <- do.call(rbind, pbmclapply(1:nsimu, function(simu_label) {
  cat(" simu", simu_label, "over", nsimu, "\r")
  network <- sample_er(p, k)
  # network <- sample_powerlaw(p)
  # network <- sample_communities(p)
  res <- Reduce("rbind", mapply(one_simu, params$n, params$covar_effect, params$covar_zi, MoreArgs = list(network = network), SIMPLIFY = FALSE))
  res$sim_label <- simu_label
  res
}, mc.cores = cores))

results_summarized <- res %>% add_column(topology = "random") %>% 
  dplyr::select(-pAUC, -AUPR) %>%  mutate(
    topology = fct_rev(factor(topology))) %>%
  mutate(method = factor(method)) %>% 
  mutate(covar_effect = fct_rev(factor(covar_effect))) %>%
  mutate(sample_size = fct_rev(factor(sample_size))) %>%
  group_by(topology, covar_effect, covar_zi, sample_size, method, sim_label) %>% 
  summarize(AUC = auc(fallout, recall), AUPR = aupr(recall, precision))

myAUC <- ggplot(results_summarized, aes(x = sample_size, y = AUC, fill = fct_reorder(method, AUC, .fun = median, .desc = TRUE))) + geom_boxplot() +
  theme_bw(base_size = 15) + scale_fill_viridis(discrete = TRUE, guide = guide_legend(title="method")) +
  theme(legend.position = c(0.75, 0.05), legend.justification = c(0.5, 0)) + xlab("") +
  facet_grid(topology ~ covar_effect, scale = "free_y")
print(myAUC)

myAUPR <- ggplot(results_summarized, aes(x = sample_size, y = AUPR, fill = fct_reorder(method, AUPR, .fun = median, .desc = TRUE))) + geom_boxplot() +
  theme_bw(base_size = 15) + scale_fill_viridis(discrete = TRUE, guide = guide_legend(title="method")) +
  theme(legend.position = c(0.1, 0.85), legend.justification = c(0.5, 0)) + xlab("") +
  facet_grid(topology ~ covar_effect, scale = "free_y")
print(myAUPR)
