library(glasso)
# library(dplyr)
library(Matrix)
library(RNAseqNet)
library(SpiecEasi)
library(PLNmodels)

sparseLLM_network <- function(Y, X) {
  res <- GLMnetwork(Y)$path
  res <- res[length(res):1]
  res
}

graphical_lasso_network <- function(Y, X, nPenalties = 50, approx = FALSE) {
  cov <- get.residuals.covariance(Y, X)
  range.penalties <- range(abs(cov[upper.tri(cov)]))
  penalties <- 10^seq(log10(range.penalties[2]), log10(max(range.penalties[1],range.penalties[2] * 1e-3)), len = nPenalties)
  out <- glasso::glassopath(cov, rholist = penalties, penalize.diagonal = FALSE, approx = approx, trace = 0)$wi
  res <- apply(out, 3,
      function(net_ ){
        net_ <- matrix(net_, ncol(cov), ncol(cov))
        Matrix(1 * net_)
      }
    )
  res <- res[length(res):1]
  res
}

graphical_lasso_covariance <- function(cov, nPenalties = 50, approx = FALSE) {
  range.penalties <- range(abs(cov[upper.tri(cov)]))
  penalties <- 10^seq(log10(range.penalties[2]), log10(max(range.penalties[1],range.penalties[2] * 1e-3)), len = nPenalties)
  out <- glasso::glassopath(cov, rholist = penalties, penalize.diagonal = FALSE, approx = approx, trace = 0)$w
  res <- apply(out, 3,
      function(cov_ ){
        cov_ <- matrix(cov_, ncol(cov), ncol(cov))
        cov_
      }
    )
  res <- res[length(res):1]
  res
}


neighborhood_selection_network <- function(X, Y, nPenalties=50) {
  graphical_lasso_network(X, Y, nPenalties=nPenalties, approx = TRUE)
}

correlation_threshold_network <- function(cov, nPenalties=100) {
  x <- cov; diag(x) <- 0
  soft_threshold <- function(lambda) {
    st <- numeric(length=length(x))
    st[which(x > lambda)]  <- x[which(x > lambda)] - lambda
    st[which(x < -lambda)] <- x[which(x < -lambda)] + lambda
    matrix(st, nrow(x), ncol(x))
  }
  path <- lapply(seq(max(abs(x)), 0, len = nPenalties), soft_threshold)
  path
}

correlation_threshold_covariance <- function(cov, nPenalties=100) {
  x <- cov; diag(x) <- 0
  soft_threshold <- function(lambda) {
    st <- numeric(length=length(x))
    st[which(x > lambda)]  <- x[which(x > lambda)] - lambda
    st[which(x < -lambda)] <- x[which(x < -lambda)] + lambda
    res <- matrix(st, nrow(x), ncol(x))
    diag(res) <- diag(cov)
    res
  }
  path <- lapply(seq(max(abs(x)), 0, len = nPenalties), soft_threshold)
  path
}

sparsePLN_network <- function(Y, X) {
  data_PLN <- prepare_data(Y, X)
  out <- PLNnetwork(Y ~ 0 + X + offset(log(Offset)), data = data_PLN, control = PLNnetwork_param(min_ratio = 1e-2, trace = 0))
  res <- sapply(out$models, function(model) {
    net_ <- model$latent_network("support")
    Matrix(1 * net_)
  })
  res
}

sparsePLN_network_nodiag <- function(Y, X) {
  data_PLN <- prepare_data(Y, X)
  out <- PLNnetwork(Y ~ 0 + X + offset(log(Offset)), data = data_PLN, control = PLNnetwork_param(min_ratio = 1e-2, trace = 0, penalize_diagonal = FALSE))
  res <- sapply(out$models, function(model) {
    net_ <- model$latent_network("support")
    Matrix(1 * net_)
  })
  res
}

sparsePLN_covariance <- function(Y, X) {
  data_PLN <- prepare_data(Y, X)
  out <- PLNnetwork(Y ~ 0 + X + offset(log(Offset)), data = data_PLN, data = data_PLN, control = PLNnetwork_param(min_ratio = 1e-2, trace = 0))
  res <- lapply(out$models, function(model) {
    model$model_par$Sigma
  })
  res
}

sparCC_network <- function(Y, X) {
  correlation_threshold_network(sparcc(Y)$Cor)
}

sparCC_covariance <- function(Y, X) {
  correlation_threshold_covariance(sparcc(Y)$Cov)
}

spiecEasi_network <- function(Y, X) {
  spiec.easi(Y, icov.select = FALSE, nlambda = 50, verbose = FALSE)$est$path
}

spiecEasi_covariance <- function(Y, X) {
  spiec.easi(Y, icov.select = FALSE, nlambda = 50, verbose = FALSE)$est$cov
}

get.residuals.covariance <- function(Y, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y))) {
  # LMs  <- lapply(1:ncol(Y), function(j) lm.fit(X, log(1 + Y[,j]), offset =  O[,j]) )
  LMs  <- lapply(1:ncol(Y), function(j) lm.fit(X, log(1 + Y[,j]), offset =  O[,j]) )
  Sigma <- cov(do.call(cbind, lapply(LMs, residuals)))
  Sigma
}
