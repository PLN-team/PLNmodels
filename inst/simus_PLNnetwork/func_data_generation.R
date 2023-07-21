library(igraph)
library(Matrix)
library(MASS)

sample_powerlaw <- function(n, power = 1, k = 1) {
  sample_pa(n, power, out.seq = rep(k, n), directed = FALSE)
}

sample_er <- function(n, k = n) {
  sample_gnm(n, k, directed = FALSE)
}

sample_regular <- function(p, k=3) {
  sample_k_regular(p, k, directed = FALSE)
}

sample_communities <- function(n, prob = c(1/2,1/4,1/4), p_in = 0.25, p_out = 0.01) {
  pref_mat <- matrix(p_out, length(prob), length(prob))
  diag(pref_mat) <- p_in
  sample_sbm(n,
             pref.matrix = pref_mat,
             block.sizes = c(rmultinom(1, n, prob)) )
}

graph2cov <- function(G, variances = rep(1, gorder(G)), u = 0.1, v = 0.3) {
  theta <- as_adj(as.undirected(G)); diag(theta) <- 0 # symmetrize the graph
  Omega <- theta * v
  diag(Omega) <- abs(min(eigen(Omega, only.values = TRUE)$values)) + u
  Sigma <- Matrix(diag(sqrt(variances)) %*% cov2cor(solve(Omega)) %*% diag(sqrt(variances)))
  symmpart(Sigma)
}

#' Simulation of community count data
#' 
#' Count data drawn under a possibily zeo inflated multinomial model with a predéfined (Gaussian) covariance structure.
#'
#' The workflow is:
#'  - draw latent abundances-basis under a centered multivariate normal with user-defined covariance matrix
#'  - use logistic transformation of abundance basis to proportions
#'  - draw counts (with a given sum N) under a multinomial distributions with previously defined proportions.
#'
#' @param n the sample size (number of communities to simulate)
#' @param mu vector of means of the latent variable (may containt covariates effect and offset); roughly equal to mean log-abundances
#' @param Sigma covariance matrix of the latent variables/species
#' @param N a vector of sequencing depth in each sample (default fixed to 3000 in all samples)
#' @param pi matrix of probabilities controling zero inflation in each sample (default fixed to 0 in all samples)
#' @return 
#' @example
#' ## Simulation settings
#' p <- 50  # number fo communities
#' n <- 100 # sample size
#' d <- 3   # number of covariates
#' k <- p   # number of edges in the network
#' network <- sample_er(p, k) # true latent network (random structure)
#' covar_effect  <- 1
#' offset_effect <- 10  ## le larger, the less effect
#' Sigma <- graph2cov(network, v = 0.3, u = 0.1)
#' ## Covariate: One-way ANOVA with 3 modalities
#' group <- factor(sample(1:d, n, replace = TRUE))
#' X <- model.matrix(~ group + 0) 
#' mu <- X %*% matrix(runif(d*p, -1,  1) , d, p)
#' ## sequencinf depth
#' N <- rnegbin(n, mu = 1000, theta = 2); N[N == 0] <- 1000
#' Y <- rMLN(n, mu = mu, Sigma = Sigma, N = N)
rMLN <- function(n, mu = matrix(0, n, p), Sigma, N = rep(3000, n), pi = rep(0, n)) {
  p <- ncol(Sigma)
  abundancies <- mu + mvrnorm(n, rep(0, p), as.matrix(Sigma))
  logistic <- function(x) { z <- exp(x); return(z / sum(z)) }
  proportions <- t(apply(abundancies, 1, logistic))
  ## sample counts
  ## first, convert matrix to list, for use with mapply
  prop.list <- lapply(seq_len(nrow(proportions)), function(i) proportions[i,])
  counts <- mapply(rmultinom, size = N, prob = prop.list, n = 1)
  counts <- t(counts)
  
}

#' @param n the sample size
#' @param mu vector of means of the latent variable (may containt covariates effect and offset)
#' @param Sigma covariance matrix of the latent variable
#' @param
rPLN <- function(n, mu = matrix(0, n, ncol(Sigma)), Sigma) {
  p <- ncol(Sigma)
  Z <- mu + mvrnorm(n, rep(0,ncol(Sigma)), as.matrix(Sigma))
  Y <- matrix(rpois(n * p, as.vector(exp(Z))),n, p)
  Y
}

