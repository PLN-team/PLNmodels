library(PLNmodels)
library(aricode)
library(ClustOfVar)
library(MASS)
library(tidyverse)
library(viridis)
theme_set(theme_bw())

rPLNblock <- function(
    n = 100,
    p = 20,
    q = 3,
    mu = rep(0, p),
    Sigma = diag(1, q, q),
    alpha = rep(1/q, q),
    depths = rep(1e4, n))  {

  if (any(is.vector(mu), ncol(mu) == 1)) {
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  }

  if (length(depths) != n) {
    depths <- rep(depths[1], n)
  }

  ## Group
  cond <- FALSE
  while(!cond) { ##
    membership <- sample.int(q, size = p, replace = TRUE, prob = alpha)
    if (length(unique(membership)) == q) cond <- TRUE
  }
  blocks <- PLNmodels:::as_indicator(membership)

  ## adjust depths
  exp_depths <- rowSums(exp(rep(1, n) %o% diag(Sigma)[membership]/2 + mu)) ## sample-wise expected depths
  offsets <- matrix(log(depths %o% rep(1, p)) - log(exp_depths), n, p)
  Z <- mu + mvrnorm(n, rep(0, q), as.matrix(Sigma)) %*% t(blocks) + offsets
  Y <- matrix(rpois(n * p, as.vector(exp(Z))), n, p)
  dimnames(Y) <- list(paste0("S", 1:n), paste0("Y", 1:p))
  rownames(offsets) <- rownames(Y)

  list(Y = Y, O = offsets, membership = membership, blocks = blocks, Sigma = Sigma)
}

### SIMULATED DATA
n <- 200
p <- 40
q <- 5
d <- 1
# R2 blanace the part of variance due to XB o Sigma
# High R2 -> part of variance due to B is important comapre to Sigma so group are harder to find
R2_target <- 0.95
params <- PLNmodels:::create_parameters(n = n, p = p, q = q, d = d, depths = 1e3)
Sigma <- toeplitz(0.75^(1:q - 1))
##  Adjusting Sigma to get a specified R2/SNR
SNR_hat <- sum(diag(cov(params$X %*% params$B))) / sum(diag(Sigma))
SNR_target <- R2_target / (1 - R2_target)
sigma2 <- SNR_hat / SNR_target
Sigma <- sigma2 * Sigma

## total number of samples should exceed max vec_n if some a drops due to lacks of count
vec_n <- c(25, 50, 75, 100)
one_simu <- function(i) {
  cat(i)
  sim <- rPLNblock(n, p, q, mu = params$X %*% params$B, Sigma = Sigma, depths = params$depths)
  data <- prepare_data(sim$Y, params$X, offset = sim$O)
  logO <- data$Offset; data$Offset <- NULL

  err_PLNblock <- sapply(vec_n,  function(n_) {
    myPLNblock <- PLNblock(Abundance ~ 0 + . + offset(logO), data = data, , subset = 1:n_, nb_blocks = q, control = PLNblock_param(trace = 0))$models[[1]]
    c(rmse_B = sqrt(mean((myPLNblock$model_par$B - params$B)^2)),
      rmse_Sigma = sqrt(mean((myPLNblock$model_par$Sigma - Sigma)^2)),
      ari = ARI(myPLNblock$membership, sim$membership)
    )
  })

  err_baseline <- sapply(vec_n,  function(n_) {
    ## revoie l'initialisation (lm  + log trasnformation)
    LMlog <- PLN(Abundance ~ 0 + . + offset(logO), data = data, , subset = 1:n_,
                 control = PLN_param(trace = 0, config_optim = list(maxeval = 1)))
    cl <- kmeans(t(LMlog$var_par$M), q)$cl
    blocks <- PLNmodels:::as_indicator(cl)
    c(rmse_B = sqrt(mean((LMlog$model_par$B - params$B)^2)),
      rmse_Sigma = sqrt(mean((LMlog$model_par$Sigma - blocks %*% Sigma %*% t(blocks))^2)),
      ari = ARI(cl, sim$membership)
    )
  })

  err_PLN <- sapply(vec_n,  function(n_) {
    myPLN <- PLN(Abundance ~ 0 + . + offset(logO), data = data, , subset = 1:n_, control = PLN_param(trace = 0))
    cl <- kmeans(t(myPLN$var_par$M), q)$cl
    blocks <- PLNmodels:::as_indicator(cl)
    c(rmse_B = sqrt(mean((myPLN$model_par$B - params$B)^2)),
      rmse_Sigma = sqrt(mean((myPLN$model_par$Sigma - blocks %*% Sigma %*% t(blocks))^2)),
      ari = ARI(cl, sim$membership)
    )
  })
  data.frame(rbind(t(err_PLNblock), t(err_PLN), t(err_baseline)),
             method = rep(c("PLNblock","PLN + kmeans", "LMlog+kmeans"), each = length(vec_n)),
             n = rep(vec_n, 3), simu = rep(i, 3*length(vec_n)))
}

res <- do.call(rbind, lapply(1:100, one_simu))

p_B <- ggplot(res) + aes(x = factor(n), y = rmse_B, fill = method) + geom_boxplot() + ylim(c(0,0.5)) +
  scale_fill_viridis_d() + ggtitle(paste("RMSE Beta (R2 =", R2_target, " p =", p, "q =",q,")"))
p_B

p_Sigma <- ggplot(res) + aes(x = factor(n), y = rmse_Sigma, fill = method) + geom_boxplot()  +
  scale_fill_viridis_d() + ylim(c(0,0.5)) +
  ggtitle(paste("RMSE Sigma (R2 =", R2_target, " p =", p, "q =",q,")"))
p_Sigma

p_ARI <- ggplot(res) + aes(x = factor(n), y = ari, fill = method) + geom_point(alpha=0.8) +
  geom_boxplot() + ylim(c(0.25,1))  + scale_fill_viridis_d() +
  ggtitle(paste("ARI (R2 =", R2_target, " p =", p, "q =",q,")"))
p_ARI

