library(PLNmodels)
data("trichoptera")
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
Y <- as.matrix(trichoptera$Abundance)
X <- model.matrix(Abundance ~ 1, data = trichoptera)
n <- nrow(Y)
p <- ncol(Y)
d <- ncol(X) # number of covariates
q <- 5 # number of PCA components
O <- matrix(0, nrow = n, ncol = p)
data <- list(Y = Y,
             X = X,
             O = O,
             w = rep(1,n))

params <- list(B = matrix(0, d, p),
               C = matrix(0, p, q),
               M = matrix(0, n, q),
               S = matrix(0.1, n, q)
        )
config <- PLNPCA_param()$config_optim

out <- PLNmodels:::nlopt_optimize_rank(data, params, config)

out$B
