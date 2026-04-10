library(tidyverse)
library(PLNmodels)
library(glmmTMB)
library(corrplot)

## 10 samples of 5 species with equal abundances, no covariance and target depths of 10,000
p <- 10
n <- 100
mu <- runif(p, -1, 1)
sd <- sqrt(sample.int(5, p, replace = TRUE))
Sigma <- diag(sd) %*% toeplitz(0.5^(1:p - 1)) %*% diag(sd)
Y <- rPLN(n, mu = mu, Sigma = Sigma, depths = c(1e3))
Y_bin <-  1*(Y > 0)

abundance_data <- Y %>% data.frame() %>% as_tibble(rownames = "ID") %>%
  pivot_longer(-ID, names_to = "var", values_to = "count", )

presence_data <- Y_bin %>% data.frame() %>% as_tibble(rownames = "ID") %>%
  pivot_longer(-ID, names_to = "var", values_to = "count")

var_names <- unique(abundance_data$var)

out_pois <- glmmTMB(count ~ 1 + (var + 0 | ID), family = poisson, data = abundance_data)

Sigma_pois <- as.matrix(VarCorr(out_pois)$cond$ID)
o <- match(var_names, gsub("var", "", rownames(Sigma_pois)))
Sigma_pois <- Sigma_pois[o,o]

out_pln <- PLN(Y ~ 1, data = abundance_data)
Sigma_pln <- sigma(out_pln)

out_bin <- glmmTMB(count ~ 1 + (var + 0 | ID), family = poisson, data = presence_data)

Sigma_bin <- as.matrix(VarCorr(out_bin)$cond$ID)
o <- match(var_names, gsub("var", "", rownames(Sigma_bin)))
Sigma_bin <- Sigma_bin[o,o]

corrplot(Sigma_pois, is.corr = FALSE, method = "color")

corrplot(Sigma_bin, is.corr = FALSE, method = "color")

corrplot(Sigma_pln, is.corr = FALSE, method = "color")

sqrt(mean((Sigma - Sigma_pln)^2))
sqrt(mean((Sigma - Sigma_pois)^2))
sqrt(mean((Sigma - Sigma_bin)^2))

