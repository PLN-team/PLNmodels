library(tidyverse)
library(PLNmodels)

create_parameters <- function(
    n = 200,
    p = 50,
    d = 2,
    rho = 0.2,
    # snr = 3,
    sigma = 1,
    depths = 100000,
    ...
) {
  ## Sigma chosen to achieve a given snr
  ## sigma <- sqrt(9*d/snr)
  list(n      = n,
       p      = p,
       X      = matrix(rnorm(n*d), nrow = n, ncol = d,
                       dimnames = list(paste0("S", 1:n), paste0("Var_", 1:d))),
       # X      = matrix(c(rep(1, n), rnorm(n*(d - 1))),
       #                 nrow = n, ncol = d,
       #                 dimnames = list(paste0("S", 1:n), paste0("Var_", 1:d))),
       Theta  = matrix(rnorm(n = p*d, sd = 1/sqrt(d)), nrow = p, ncol = d),
       Sigma  = sigma * toeplitz(x = rho^seq(0, p-1)),
       depths = depths)
}

params <- create_parameters()
Theta <- params$Theta

## Extract X
X <- params$X

## Extract Y
Y <- rPLN(n = nrow(X), mu = tcrossprod(X, Theta), Sigma = params$Sigma, depths = params$depths)

data <- prepare_data(Y, X, offset = "none")
O <- rowSums(Y)
model <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data,
             control = list(trace = 0, covariance = "fixed", prec_matrix = solve(params$Sigma)))
Theta_hat <- coef(model)
model$get_vcov_hat("wald", Y, X)
Theta_se_wald <- standard_error(model)
model$get_vcov_hat("louis", Y, X)
Theta_se_louis <- standard_error(model)
model$get_vcov_hat("sandwich", Y, X)
Theta_se_sandwich <- standard_error(model)

data.frame(
  Theta = rep(c(Theta), 3),
  Theta_hat = rep(c(Theta_hat), 3),
  se = c(Theta_se_wald, Theta_se_louis, Theta_se_sandwich),
  method = rep(c("wald", "louis", "sandwich"), each = length(c(Theta))) ) %>%
  ggplot(aes(x = Theta, y = Theta_hat)) +
  geom_errorbar(aes(ymin = Theta_hat - 2 * se,
                    ymax = Theta_hat + 2 * se), color = "blue") + facet_wrap(~ method) +
  geom_abline(slope = 1, intercept = 0) + labs(x = "True value", y = "Mean estimate") + theme_bw() -> p

print(p)
