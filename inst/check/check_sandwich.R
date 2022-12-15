library(tidyverse)
library(PLNmodels)

nb_cores <- 10
options(future.fork.enable = TRUE)

params <- PLNmodels:::create_parameters()
Theta <- params$Theta

## Extract X
X <- params$X

## Extract Y
Y <- rPLN(n = nrow(X), mu = tcrossprod(X, Theta), Sigma = params$Sigma, depths = params$depths)

data <- prepare_data(Y, X, offset = "none")
O <- rowSums(Y)


future::plan("multicore", workers = nb_cores)
model <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data, control = PLN_param(jackknife = TRUE, bootstrap = TRUE))
future::plan("sequential")


Theta_hat <- coef(model)
Theta_se_var <- standard_error(model)
Theta_se_jk <- standard_error(model, "jackknife")
Theta_se_bt <- standard_error(model, "bootstrap")

data.frame(
  Theta = rep(c(Theta), 3),
  Theta_hat = rep(c(Theta_hat), 3),
  se = c(Theta_se_var, Theta_se_jk, Theta_se_bt),
  method = rep(c("variational", "jackknife", "bootstrap"), each = length(c(Theta))) ) %>%
  ggplot(aes(x = Theta, y = Theta_hat)) +
  geom_errorbar(aes(ymin = Theta_hat - 2 * se,
                    ymax = Theta_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)

