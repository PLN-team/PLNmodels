library("PLNmodels")

# choose data
n <- 100
p <- 2
counts <- matrix(rpois(n*p, c(5,11)), n, p)
covariates <- matrix(1, n, 1)

# Fit PLN on data
data  <- prepare_data(counts, covariates)
model <- PLN(Abundance ~ 1 + offset(log(Offset)), data = data)

# take training data
# get parameters for the VE-step
new_n <- n
ind <- 1:new_n
new_data <- data[ind, , drop=FALSE]
new_responses  <- new_data$Abundance
new_covariates <- matrix(new_data$V1, ncol=1)
new_offsets    <- matrix(rep(log(new_data$Offset), p), ncol = p)
new_weights <- rep(1, new_n)

B <- model$model_par$B
Omega <- model$model_par$Omega
M_init <- model$var_par$M[ind, , drop=FALSE]
S_init <- model$var_par$S[ind, , drop=FALSE]

M_init <- matrix(0  , new_n, p)
S_init <- matrix(0.1, new_n, p)

args <- list(data = list(Y = new_responses, X = new_covariates, O = new_offsets, w = new_weights),
             ## Initialize the variational parameters with the new dimension of the data
             params = list(M = M_init, S = S_init),
             B = as.matrix(B),
             Omega = as.matrix(Omega) ,
             config = PLN_param()$config_optim)
VE <- do.call(PLNmodels:::nlopt_optimize_vestep, args)

mse <- function(a, b) sum((a-b)^2)
print(mse(VE$S**2, S_init**2))
print(mse(VE$M, M_init))
print(mse(VE$Ji, model$loglik_vec[ind]))
