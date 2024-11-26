library(tidyverse)
library(PLNmodels)
library(future.apply)
set.seed(1234)

plan(multisession, workers = 10)

rmse <- function(theta_hat, theta_star) {
  sqrt(mean((theta_hat - theta_star)^2))
}

params <- PLNmodels:::create_parameters(n = 250, p = 10, d = 1, depths = 1e3)
X <- params$X
B <- params$B
conf <- list(variational_var = TRUE, jackknife = FALSE, bootstrap = 50, sandwich_var = TRUE)
seq_n <- c(50, 100, 250)

one_simu <- function(s) {
  cat("\nsimu", s)
  do.call(rbind, lapply(seq_n, function(n) {
    X_sub <- X[1:n, , drop = FALSE]
    Y <- rPLN(n = n, mu = X_sub %*% B, Sigma = params$Sigma, depths = params$depths)
    data <- prepare_data(Y, X_sub, offset = "none")
    logO <- attr(Y, "offsets")
    model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(trace = FALSE, config_post = conf))

    B_hat <- coef(model)
    vcov_sandwich    <- attr(coef(model), "vcov_sandwich")
    vcov_bootstrap   <- attr(coef(model), "vcov_bootstrap")
    vcov_variational <- attr(coef(model), "vcov_variational")

    data.frame(rmse = rmse(B_hat, B),
               cover_sandwich    = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_sandwich))) < 1.96),
               cover_bootstrap   = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_bootstrap))) < 1.96),
               cover_variational = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_variational))) < 1.96),
               sample_size = n,
               simu = s)
  }))
}

res <- do.call(rbind, lapply(1:50, one_simu))

res %>%
  pivot_longer(-c(rmse, sample_size, simu), values_to = "coverage", names_to = "estimator") %>%
  ggplot() + aes(x = estimator, y = coverage, fill = factor(sample_size)) + geom_boxplot()

### Single test
Y <- rPLN(n = nrow(X), mu = X %*% B, Sigma = params$Sigma, depths = params$depths)
data <- prepare_data(Y, X, offset = "none")
logO <- attr(Y, "offsets")
model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(config_post = conf))

B_hat <- coef(model)
B_se_var <- standard_error(model, "variational")
B_se_bt  <- standard_error(model, "bootstrap")
B_se_sw  <- standard_error(model, "sandwich")

data.frame(
  B = rep(c(B), 3),
  B_hat = rep(c(B_hat), 3),
  se = c(B_se_var, B_se_bt, B_se_sw),
  method = rep(c("variational", "bootstrap", "sandwich"), each = length(c(B))) ) %>%
  ggplot(aes(x = B, y = B_hat)) +
  geom_errorbar(aes(ymin = B_hat - 2 * se,
                    ymax = B_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)

