library(tidyverse)
library(PLNmodels)
set.seed(1234)

rmse <- function(theta_hat, theta_star) {
  sqrt(mean((theta_hat - theta_star)^2))
}

params <- PLNmodels:::create_parameters(n = 100, p = 10, d = 1, depths = 1e3)
X <- params$X
B <- params$B
conf <- list(variational_var = TRUE, jackknife = TRUE, bootstrap = FALSE, sandwich_var = TRUE)


one_simu <- function(s) {

  Y <- rPLN(n = nrow(X), mu = X %*% B, Sigma = params$Sigma, depths = params$depths)
  data <- prepare_data(Y, X, offset = "none")
  logO <- attr(Y, "offsets")
  model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(trace = FALSE, config_post = conf))

  B_hat <- coef(model)
  vcov_sandwich    <- attr(coef(model), "vcov_sandwich")
  vcov_jackknife   <- attr(coef(model), "vcov_sandwich")
  vcov_variational <- attr(coef(model), "vcov_variational")

  data.frame(rmse = rmse(B_hat, B),
       cover_sandwich    = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_sandwich))) < 1.96),
       cover_jackknife   = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_jackknife))) < 1.96),
       cover_variational = mean(abs(as.numeric(B_hat - B) %*% solve(chol(vcov_variational))) < 1.96),
       simu = s)
}

res <- do.call(rbind, lapply(1:50, one_simu))

boxplot(res$cover_sandwich, res$cover_jackknife, res$cover_variational)

### Single test

B_se_var <- standard_error(model, "variational")
B_se_jk  <- standard_error(model, "jackknife")
B_se_sw  <- standard_error(model, "sandwich")

Y <- rPLN(n = nrow(X), mu = X %*% B, Sigma = params$Sigma, depths = params$depths)
data <- prepare_data(Y, X, offset = "none")
logO <- attr(Y, "offsets")
model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(config_post = conf))

data.frame(
  B = rep(c(B), 3),
  B_hat = rep(c(B_hat), 3),
  se = c(B_se_var, B_se_jk, B_se_sw),
  method = rep(c("variational", "jackknife", "sandwich"), each = length(c(B))) ) %>%
  ggplot(aes(x = B, y = B_hat)) +
  geom_errorbar(aes(ymin = B_hat - 2 * se,
                    ymax = B_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)

