library(tidyverse)
library(PLNmodels)
set.seed(1234)

nb_cores <- 10
options(future.fork.enable = TRUE)

params <- PLNmodels:::create_parameters(n = 50, p = 10, d = 1, depths = 1e3)
X <- params$X
B <- params$B
Y <- rPLN(n = nrow(X), mu = X %*% B, Sigma = params$Sigma, depths = params$depths)

data <- prepare_data(Y, X, offset = "none")
logO <- attr(Y, "offsets")

conf <- list(variational_var = TRUE, jackknife = TRUE, bootstrap = nrow(Y), sandwich_var = TRUE)
future::plan("multicore", workers = nb_cores)
model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(config_post = conf))
future::plan("sequential")

B_hat <- coef(model)
B_se_var <- standard_error(model, "variational")
B_se_jk  <- standard_error(model, "jackknife")
B_se_bt  <- standard_error(model, "bootstrap")
B_se_var <- standard_error(model, "variational")
B_se_jk  <- standard_error(model, "jackknife")
B_se_sw  <- standard_error(model, "sandwich")

data.frame(
  B = rep(c(B), 4),
  B_hat = rep(c(B_hat), 4),
  se = c(B_se_var, B_se_jk, B_se_bt, B_se_sw),
  method = rep(c("variational", "jackknife", "bootstrap", "sandwich"), each = length(c(B))) ) %>%
  ggplot(aes(x = B, y = B_hat)) +
  geom_errorbar(aes(ymin = B_hat - 2 * se,
                    ymax = B_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)

