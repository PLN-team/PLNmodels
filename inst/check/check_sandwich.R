library(tidyverse)
library(PLNmodels)

params <- PLNmodels:::create_parameters()
Theta <- params$Theta

## Extract X
X <- params$X

## Extract Y
Y <- rPLN(n = nrow(X), mu = tcrossprod(X, Theta), Sigma = params$Sigma, depths = params$depths)

data <- prepare_data(Y, X, offset = "none")
O <- rowSums(Y)
model <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data,
             control = PLN_param(trace = 0, covariance = "fixed", Omega = solve(params$Sigma)))

Theta_hat <- coef(model)
Theta_se_var <- standard_error(model)

model$variance_jackknife(Abundance ~ 0 + . + offset(log(O)), data = data)
Theta_se_jk <- standard_error(model, "jackknife") * sqrt(nrow(X))

# model$vcov_model("sandwich", Y, X)
# Theta_se_sandwich <- standard_error(model)

data.frame(
  Theta = rep(c(Theta), 2),
  Theta_hat = rep(c(Theta_hat), 2),
  se = c(Theta_se_var, Theta_se_jk),
  method = rep(c("variational", "jackknife"), each = length(c(Theta))) ) %>%
  ggplot(aes(x = Theta, y = Theta_hat)) +
  geom_errorbar(aes(ymin = Theta_hat - 2 * se,
                    ymax = Theta_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)
