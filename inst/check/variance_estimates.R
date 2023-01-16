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

conf <- list(variational_var = TRUE, jackknife = TRUE, bootstrap = nrow(Y))
future::plan("multicore", workers = nb_cores)
model <- PLN(Abundance ~ 0 + . + offset(logO), data = data, control = PLN_param(config_post = conf))
future::plan("sequential")

B_hat <- coef(model)
B_se_var <- standard_error(model, "variational")
B_se_jk  <- standard_error(model, "jackknife")
B_se_bt  <- standard_error(model, "bootstrap")

data.frame(
  B = rep(c(B), 3),
  B_hat = rep(c(B_hat), 3),
  se = c(B_se_var, B_se_jk, B_se_bt),
  method = rep(c("variational", "jackknife", "bootstrap"), each = length(c(B))) ) %>%
  ggplot(aes(x = B, y = B_hat)) +
  geom_errorbar(aes(ymin = B_hat - 2 * se,
                    ymax = B_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)

